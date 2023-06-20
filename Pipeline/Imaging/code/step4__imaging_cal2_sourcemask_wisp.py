from glob import glob
import os
import urllib
import shutil
from astropy.io import fits
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.stats import SigmaClip, sigma_clipped_stats
from scipy import ndimage
from astropy.coordinates import SkyCoord
from astropy import units as u
from stcal.dqflags import interpret_bit_flags
from jwst.resample.resample_utils import build_mask
from scipy import ndimage
from photutils.segmentation import detect_sources, SourceCatalog
from photutils.utils import circular_footprint
from astropy.stats import median_absolute_deviation, biweight_location, biweight_scale
from jwst.datamodels import ImageModel
import warnings
import yaml
warnings.filterwarnings("ignore")
#################################################################################
def measure_moments(imgdata, inthreshmask, kron_threshes=[5,12,20,40], ndilate=[3, 5, 10, 40]):
    threshmask = np.copy(inthreshmask)
    #objmask = Zlabeled == objid
    #moments
    Zlabeled, Nlabels = ndimage.label(threshmask)
    x = np.arange(0, imgdata.shape[1])
    y = np.arange(0, imgdata.shape[0])
    xv, yv = np.meshgrid(x, y)
    
    r_krons = np.zeros(Nlabels+1, dtype=float)
    kron_map = np.zeros(imgdata.shape, dtype=float)
    for objid in range(1,Nlabels+1):   
        #print(objid)
        objmask = Zlabeled == objid
        cent_x = np.nansum(imgdata[(objmask)] * xv[(objmask)]) / np.nansum(imgdata[(objmask)])
        cent_y = np.nansum(imgdata[(objmask)] * yv[(objmask)]) / np.nansum(imgdata[(objmask)])
        #sigmas
        #sig_x = np.nansum(imgdata[(objmask)] * ((xv-cent_x)**2)[(objmask)]) / np.nansum(imgdata[(objmask)])
        #sig_y = np.nansum(imgdata[(objmask)] * ((yv-cent_y)**2)[(objmask)]) / np.nansum(imgdata[(objmask)])
        #x2 =  np.nansum(imgdata[(objmask)] * ((xv)**2)[(objmask)]) / np.nansum(imgdata[(objmask)]) - cent_x**2
        #kron rad
        rpix = np.sqrt((xv-cent_x)**2 + (yv-cent_y)**2)
        r_krons[objid] = np.nansum((imgdata * rpix)[(objmask)]) / np.nansum(imgdata[(objmask)])
        kron_map[(objmask)] = r_krons[objid]
     
    mask_cube = np.zeros([len(kron_threshes), imgdata.shape[0], imgdata.shape[1]], dtype=bool)
    for it in range(len(kron_threshes)):
        mask_cube[it,:,:] = dilate_mask(kron_map > kron_threshes[it], niter = ndilate[it])
    return np.any(mask_cube, axis=0), mask_cube
#################################################################################
def fast_minarea(inmask, minpix=20):
    threshmask = np.copy(inmask)
    # now identify the objects and remove those above a threshold
    Zlabeled, Nlabels = ndimage.label(threshmask)
    label_bins = np.arange(0, Nlabels+2)
    label_size, label_bins = np.histogram(Zlabeled, bins=label_bins)
    #label_size = [(Zlabeled == label).sum() for label in range(Nlabels + 1)] 
    bad_labels = label_bins[:-1][(label_size < minpix)]
    for label in bad_labels: threshmask[Zlabeled == label] = False
    return threshmask
#################################################################################
def dilate_mask(threshmask, niter=1):
    #dilation
    grow1  = ndimage.generate_binary_structure(2, 1)
    grow2  = ndimage.generate_binary_structure(2, 2)
    threshmask_exp = np.copy(threshmask)
    for i in range(niter):
        threshmask_exp = ndimage.binary_dilation(threshmask_exp, structure=grow1)
        threshmask_exp = ndimage.binary_dilation(threshmask_exp, structure=grow2)
    return threshmask_exp
#################################################################################
def threshmask_dict(imgdata, inerr, median, std, dq_mask=None, layerlist=[dict(sigma=5.0, thresh=3.0, minpix=1000, ndilate='kron'),
                                                      dict(sigma=3.0, thresh=3.0, minpix=None, ndilate=3),
                                                      dict(sigma=1.0, thresh=2.0, minpix=20,   ndilate=1)]):
    imgcut = np.copy(imgdata)
    #fix zero errors in err
    err = np.copy(inerr)
    err[(err == 0.0)] = np.nanmedian(err)
    #layerlist: list of dicts of mask layer e.g. [dict(sigma=5.0, thresh=3.0, minpix=20, ndilate=10)]
    #measure stats
    #mean,  median,  std = sigma_clipped_stats(imgcut)
    #loop around layers
    mask_array = np.zeros([len(layerlist), imgcut.shape[0], imgcut.shape[1]], dtype=bool)
    for il, layer in enumerate(layerlist):

        #make big tier
        imgcut_smth = ndimage.gaussian_filter(imgcut, layer['sigma'])
        sigma_smth = (imgcut_smth - median)/err
        thresh_mask = sigma_smth > layer['thresh']
        #minpix
        if layer['minpix'] is not None:
            thresh_mask = fast_minarea(thresh_mask, minpix=layer['minpix'])
            
        if layer['ndilate'] == 'kron':
            threshmask_exp, mask_cube = measure_moments(imgcut, thresh_mask)
        else:
            #dilate
            threshmask_exp = dilate_mask(thresh_mask, niter=layer['ndilate'])
        #store
        mask_array[il,:,:] = threshmask_exp

    #stack
    return np.any(mask_array, axis=0), mask_array
#################################################################################
def fast_mosaic_mask(mosaic, FILTER, save_dir):
    #imheader = cal_hdu[1].header
    wcs_mosaic =WCS(mosaic['SCI'].header)
    mosaic_mask = np.zeros(mosaic['SCI'].data.shape, dtype=bool)
    #detection thresh scale
    mean,  median,  std = sigma_clipped_stats(mosaic['SCI'].data, mask=(mosaic['SCI'].data==0))
    #wcs=WCS(imheader) #see above
    for ii in range(np.ceil(mosaic['SCI'].data.shape[0]/2000).astype(int)):
        for jj in range(np.ceil(mosaic['SCI'].data.shape[1]/2000).astype(int)):
    
            #mask_data=mosaic['SCI'].data
            #ra_pix, dec_pix = wcs_mosaic.all_pix2world([0,0,2048,2048],[0,2048,2048,0],0) ##Changed 2048 to 2038 to not mask borders
            #x_mos,y_mos = wcs_mosaic.all_world2pix(ra_pix,dec_pix,0)
            
            keep_area = [[np.max([0,ii*2000]), np.min([mosaic['SCI'].data.shape[0],(ii+1)*2000])],
                        [np.max([0,jj*2000]), np.min([mosaic['SCI'].data.shape[1],(jj+1)*2000])]]    
            cutout_area = [[np.max([0,keep_area[0][0]-500]), np.min([mosaic['SCI'].data.shape[0],keep_area[0][1]+500])],
                            [np.max([0,keep_area[1][0]-500]), np.min([mosaic['SCI'].data.shape[1],keep_area[1][1]+500])]]

            mosaic_cutout = mosaic['SCI'].data[cutout_area[0][0]:cutout_area[0][1], cutout_area[1][0]:cutout_area[1][1]]
            err_cutout   = mosaic['ERR'].data[cutout_area[0][0]:cutout_area[0][1], cutout_area[1][0]:cutout_area[1][1]]

            view = slice(cutout_area[0][0],cutout_area[0][1] ), slice(cutout_area[1][0],cutout_area[1][1])		#check order
            wcs_cut = wcs_mosaic.slice(view)

            #make mask 
            mask_cutout, mask_cube = threshmask_dict(mosaic_cutout, err_cutout, median, std)
            #get corners and mosaic cutout
            mosaic_mask[keep_area[0][0]:keep_area[0][1], keep_area[1][0]:keep_area[1][1]] = mask_cutout[
            keep_area[0][0]-cutout_area[0][0]:(keep_area[0][1] - keep_area[0][0])+(keep_area[0][0]-cutout_area[0][0]), 
            keep_area[1][0]-cutout_area[1][0]:(keep_area[1][1] - keep_area[1][0])+(keep_area[1][0]-cutout_area[1][0])] 
            print(ii, jj)

    mask_hdu = fits.PrimaryHDU(mosaic_mask.astype(float), header=mosaic['SCI'].header, do_not_scale_image_data=True)
    mask_hdu.writeto(save_dir+'sourcemask_cal2_v1.fits')
    return mosaic_mask, mask_hdu
#################################################################################
def make_src_mask(cal_hdu, mask_hdu):
    imheader = cal_hdu[1].header
    wcs=WCS(imheader) #see above

    wcs_mask=WCS(mask_hdu['SCI'].header)
    mask_data=mask_hdu['SCI'].data
    mask=np.zeros((2048,2048))
    for j in range(2048):
            ra_pix,dec_pix=wcs.all_pix2world(np.arange(2048),np.zeros(2048)+j,0) ##Changed 2048 to 2038 to not mask borders
            x_mask,y_mask=wcs_mask.all_world2pix(ra_pix,dec_pix,0)
            if np.any((x_mask > mask_data.shape[1])|(y_mask > mask_data.shape[0])):
                print(np.max(x_mask), np.max(y_mask), mask_data.shape, j)
            y_round = np.array(np.round(y_mask),dtype='int')
            y_round[(y_round >= mask_data.shape[0])] = mask_data.shape[0]-1
            x_round = np.array(np.round(x_mask),dtype='int')
            x_round[(x_round >= mask_data.shape[1])] = mask_data.shape[1]-1
            mask[j,:]=mask_data[y_round,x_round]

    return mask == 1
#################################################################################  
def make_wisp_template(FILTER, data_dir, save_dir, mask_hdu):
    if FILTER in ['F277W', 'F356W', 'F410M', 'F444W']: 
        camlist = ['long']
    else:
        camlist = [1,2,3,4]

    for module in ['a','b']:
        for camera in camlist:
            fits_files_cal = glob(data_dir+'jw*nrc%s%s_crf.fits'%(module,camera))
            all_cal  = np.zeros([len(fits_files_cal),2048,2048], dtype=float)
            all_mask = np.zeros([len(fits_files_cal),2048,2048], dtype=bool)
            for i, cal_filename in enumerate(fits_files_cal):
                cal_hdu = fits.open(cal_filename)
                all_cal[i,:,:] = cal_hdu['SCI'].data
                all_mask[i,:,:] = make_src_mask(cal_hdu, mask_hdu)
                #med sub
                #mean, med, stddev = sigma_clipped_stats(all_cal[i,:,:], mask=all_mask[i,:,:], sigma=3.0, axis=0)
                all_cal[i,:,:] = all_cal[i,:,:] 

            mask_mean, mask_median, stddev = sigma_clipped_stats(all_cal, mask=all_mask, sigma=3.0, axis=0)
            fits.writeto(save_dir+'cal2_srcmask_median_stack_nrc%s%s.fits'%(module,camera), mask_median, overwrite=True)
#################################################################################
def main():
    #load reduction parameters
    if not os.path.isfile('reduction_params.yml'):
        sys.exit('no reduction_params.yml file found')

    with open('reduction_params.yml') as f:
    # use safe_load instead load
        dataMap = yaml.safe_load(f)

    os.environ["CRDS_CONTEXT"] = dataMap['pmap']
    basedir                    = dataMap['basedir']

    FILTERS = ['F115W', 'F200W', 'F356W']
    for FILTER in FILTERS:
        save_dir = basedir+FILTER+'/mycals/'
        pipe3_dir = basedir+FILTER+'/pipe3_skyfix/'
        pipe2_dir = basedir+FILTER+'/pipe2_mywcs/'

        #load mosaic
        mosaic_name =  "%sstack_%s_pipe3.fits"%(pipe3_dir,FILTER)
        mosaic = fits.open(mosaic_name)

        #make mask
        #mosaic_mask = fast_mosaic_mask(mosaic, FILTER, save_dir)
        mask_hdu = fits.open(save_dir+'sourcemask_cal2_v1.fits')

        #remake wisps
        make_wisp_template(FILTER, pipe2_dir, save_dir, mask_hdu)
##########################################################  
if __name__ == "__main__":
    main()
