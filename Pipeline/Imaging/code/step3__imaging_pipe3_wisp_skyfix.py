import matplotlib.pyplot as plt
from scipy.stats import binned_statistic
import os
from jwst.pipeline import Image3Pipeline
from astropy.visualization import simple_norm
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table, Column, join
from matplotlib.gridspec import GridSpec
#from regions import Regions
from scipy.stats import gaussian_kde
from matplotlib.colors import Normalize 
#from scipy.interpolate import interp
from matplotlib import cm
from astropy.stats import SigmaClip, sigma_clipped_stats
from scipy import interpolate
from glob import glob
from astropy.visualization import ZScaleInterval
from scipy import ndimage
from photutils.segmentation import detect_sources, SourceCatalog
from photutils.utils import circular_footprint
from astropy.stats import median_absolute_deviation, biweight_location, biweight_scale
from jwst.resample.resample_utils import build_mask
from multiprocessing import Pool
from jwst.datamodels import ImageModel
from datetime import datetime
import yaml
#################################################################################
def contiguity(segment_map, iobj):
    obj_mask = (segment_map.data == iobj)
    osize =  obj_mask.shape
    n_connected = np.zeros(obj_mask.shape, dtype=int)
    n_connected[0+0:osize[0]-1,0+0:osize[1]-0] += obj_mask[0+1:osize[0]-0,0+0:osize[1]-0].astype(int)
    n_connected[0+1:osize[0]-0,0+0:osize[1]-0] += obj_mask[0+0:osize[0]-1,0+0:osize[1]-0].astype(int)
    n_connected[0+0:osize[0]-0,0+0:osize[1]-1] += obj_mask[0+0:osize[0]-0,0+1:osize[1]-0].astype(int)
    n_connected[0+0:osize[0]-0,0+1:osize[1]+0] += obj_mask[0+0:osize[0]-0,0+0:osize[1]-1].astype(int)
    #sum flux in obj
    contig_frac = np.mean((n_connected[(obj_mask)]).astype(float))
    return contig_frac
#################################################################################
def make_snowball_mask(dq_data, lzero_cont_skip=True):

    dq_mask = (build_mask(dq_data,  '~JUMP_DET') == 0)
    segment_map = detect_sources(dq_mask, 0.5, npixels=50)
    cat = SourceCatalog(dq_mask, segment_map)
    tab =  cat.to_table()
    
    layers = [[50, 0.5, 5], [200, 1.0, 10], [500, 1.0, 20], [1000, 1.0, 30]]
    save_mask =np.zeros([len(layers), dq_mask.shape[0], dq_mask.shape[1]], dtype=bool)
    for il, layer in enumerate(layers):
        seg = segment_map.copy()
        cal_sel = tab[(tab['eccentricity'].value < layer[1])&(tab['area'].value > layer[0])]

        #new step which checks large objects are contigious
        if (il == 0)&(lzero_cont_skip):
            cal_sel2 = cal_sel
        else:
            cont_rate = np.zeros(len(cal_sel), dtype=float)
            for io, obj in enumerate(cal_sel): cont_rate[io] =  contiguity(segment_map, obj['label'])
            cal_sel2 = cal_sel[(cont_rate > 3)]

        seg.keep_labels(labels=cal_sel2['label'].data)
        footprint = circular_footprint(radius=layer[2])
        save_mask[il,:,:] = seg.make_source_mask(footprint=footprint)

    return np.any(save_mask, axis=0)
#################################################################################
def photutils_mask(imgdata, thresh=2.0, npix=10, sigma=2.0, radius=3):
    mean, median, stddev = sigma_clipped_stats(imgdata, sigma=3.0)
    imgcut = np.copy(imgdata)
    imgcut_smth = ndimage.gaussian_filter(imgcut, sigma)
    segment_map = detect_sources(imgcut_smth, thresh*stddev, npixels=npix)
    footprint = circular_footprint(radius=radius)
    save_mask = segment_map.make_source_mask(footprint=footprint)

    #mask for bigger objs
    big_segment_map = detect_sources(imgcut_smth, thresh*stddev, npixels=500)
    if big_segment_map == None: 
        return save_mask
    big_footprint = circular_footprint(radius=10)
    big_save_mask = big_segment_map.make_source_mask(footprint=big_footprint)
    return (save_mask)|(big_save_mask)
#################################################################################
def proc_single_cal(params):
    filename      = params['filename']
    wisp_template = params['wisp_template']
    save_dir      = params['save_dir']

    filtname = save_dir + (filename.split('/')[-1])
    #cal_hdu = fits.open(cal_filename)
    dm = ImageModel(filename)
    original_image = np.copy(dm.data)

    #apply wisp template
    #if FILTER in ['F115W','F150W','F200W']:
    dm.data -= wisp_template

    #make sball mask
    sball_mask = make_snowball_mask(dm.dq)
    dm.dq[(sball_mask)] = 1

    #set median sky to zero
    dq_mask = (build_mask(dm.dq,  '~DO_NOT_USE+NON_SCIENCE')==0)
    mean, median, stddev = sigma_clipped_stats(dm.data, mask=dq_mask, sigma=3.0)
    dm.data -= median
    #update meta data
    dm.meta.background.level = 0.0
    dm.meta.background.subtracted = True

    dm.data[(original_image == 0.)] = 0.

    dm.write(filtname, overwrite=True) #.replace('cal.fits','cal.fits')
#################################################################################
def cal_wisp_sub_para(save_dir, data_dir, cal_dir, FILTER, n_procs=40):
    if FILTER in ['F277W', 'F356W', 'F410M', 'F444W']: 
        camlist = ['long']
    else:
        camlist = [1,2,3,4]
    for module in ['a','b']:
        for camera in camlist:
            if camera != 'long':
                wisp_template_raw = fits.getdata(cal_dir + 'pipe3_masked_median_stack_nrc%s%s.fits'%(module,camera))
                biwt = biweight_location(wisp_template_raw, ignore_nan=True)
                wisp_template = wisp_template_raw - biwt
            else:
                wisp_template = np.zeros([2048,2048], dtype=float)

            fits_files_cal = glob(data_dir+'jw*nrc%s%s_crf.fits'%(module,camera))

            with Pool(n_procs) as pool:
                dictlist = [dict(filename=fname, wisp_template=wisp_template, save_dir=save_dir) for fname in fits_files_cal]
                pool.map(proc_single_cal, dictlist)
#################################################################################
def restack(params):
    #load wcs 
    FILTER = params['FILTER']
    basedir = params['basedir']

    wcs_pickle = basedir + 'output_grid_wcs.npz'
    save_wcs = np.load(wcs_pickle)

    #Step 3:
    save_dir = basedir+FILTER+'/pipe3_skyfix/'
    fits_files_corr=glob(save_dir + 'jw*crf.fits')

    pipe = Image3Pipeline()
    pipe.output_dir=save_dir
    pipe.tweakreg.skip=True
    #pipe.tweakreg.align_to_gaia=False
    pipe.outlier_detection.skip=True
    pipe.save_results=True
    pipe.skymatch.skip=True
    #fixed grid for resample
    pipe.resample.pixel_scale = 0.03
    pipe.resample.output_shape = save_wcs['output_shape'].tolist()           #nx first and ny
    pipe.resample.crpix = save_wcs['crpix'].tolist()
    pipe.resample.crval = save_wcs['crval'].tolist()
    pipe.resample.rotation = save_wcs['rotation'].tolist()
    pipe.run(fits_files_corr)

    os.system('mv %sstep_i2d.fits %sstack_%s_pipe3.fits'%(save_dir, save_dir, FILTER))
#################################################################################
def main():
    #load reduction parameters
    if not os.path.isfile('reduction_params.yml'):
        sys.exit('no reduction_params.yml file found')

    with open('reduction_params.yml') as f:
    # use safe_load instead load
        dataMap = yaml.safe_load(f)

    os.environ["CRDS_SERVER_URL"]   = "https://jwst-crds.stsci.edu"
    os.environ["CRDS_DATA"]         = dataMap['cache']
    os.environ["CRDS_PATH"]         = dataMap['cache']
    os.environ["CRDS_CONTEXT"]      = dataMap['pmap']
    basedir                         = dataMap['basedir']

    #first loop over filters and make wisp templates
    FILTERS = ['F115W', 'F200W', 'F356W']
    for FILTER in FILTERS:
        data_dir = basedir+FILTER+'/pipe2_mywcs/'
        save_dir = basedir+FILTER+'/pipe3_skyfix/'
        cal_dir = basedir+FILTER+'/mycals/'
        os.chdir(save_dir)

        if FILTER in ['F277W', 'F356W', 'F410M', 'F444W']: 
            camlist = ['long']
        else:
            camlist = [1,2,3,4]

        for module in ['a','b']:
            for camera in camlist:
                fits_files_cal = glob(data_dir+'jw*nrc%s%s_crf.fits'%(module,camera))
                #load images
                all_cal = np.array([fits.getdata(cal_filename) for cal_filename in fits_files_cal])
                #sigma clipped stack
                simp_mean, simp_median, simp_stddev = sigma_clipped_stats(all_cal, sigma=3.0, axis=0)


                all_mask = np.zeros(all_cal.shape, dtype=float)
                for i in range(all_cal.shape[0]):
                    my_fix = all_cal[i,:,:] - simp_median
                    my_fix_sub = my_fix - np.nanmean(my_fix)        #should be meadian
                    all_mask[i,:,:] = photutils_mask(my_fix_sub)

                mask_mean, mask_median, stddev = sigma_clipped_stats(all_cal, mask=all_mask, sigma=3.0, axis=0)

                mask_median[(np.all(all_mask, axis=0))] = simp_median[(np.all(all_mask, axis=0))]
                mask_mean[(np.all(all_mask, axis=0))] = simp_median[(np.all(all_mask, axis=0))]

                fits.writeto(cal_dir+'pipe3_median_stack_nrc%s%s.fits'%(module,camera),simp_median,overwrite=True)
                fits.writeto(cal_dir+'pipe3_masked_median_stack_nrc%s%s.fits'%(module,camera),simp_median,overwrite=True)

    #now file all and sky sub
    for FILTER in FILTERS:
        #now fix crf cal files     
        data_dir = basedir+FILTER+'/pipe2_mywcs/'
        save_dir = basedir+FILTER+'/pipe3_skyfix/'
        cal_dir = basedir+FILTER+'/mycals/'
        cal_wisp_sub_para(save_dir, data_dir, cal_dir, FILTER)

    # #now stack
    # with Pool(3) as pool:
    #     dictlist = [dict(FILTER=FILT, basedir=basedir) for FILT in FILTERS]
    #     pool.map(restack, dictlist)

    restack(dict(FILTER='F356W', basedir=basedir))
    restack(dict(FILTER='F200W', basedir=basedir))
    restack(dict(FILTER='F115W', basedir=basedir))
##########################################################  
if __name__ == "__main__":
    main()
