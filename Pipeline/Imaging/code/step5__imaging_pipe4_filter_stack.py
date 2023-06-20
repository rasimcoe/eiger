from glob import glob
import os
import urllib
os.environ["CRDS_DATA"] = "/scratch/mruari/EIGER/cache/crds_cache"
os.environ["CRDS_PATH"] = "/scratch/mruari/EIGER/cache/crds_cache"
os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"
import shutil
from astropy.io import fits
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.stats import SigmaClip, sigma_clipped_stats, biweight_location
from scipy import ndimage
from jwst.pipeline import Image3Pipeline
from astropy.coordinates import SkyCoord
from astropy import units as u
from stcal.dqflags import interpret_bit_flags
from jwst.resample.resample_utils import build_mask
from scipy import ndimage
from photutils.segmentation import detect_sources, SourceCatalog
from photutils.utils import circular_footprint
from jwst.datamodels import ImageModel
from datetime import datetime
import yaml
from multiprocessing import Pool
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
def make_src_mask_dm(cal_dm, mask_data, mask_hdr):
    #imheader = cal_hdu[1].header
    wcs = cal_dm.meta.wcs #see above

    wcs_mask=WCS(mask_hdr)
    #mask_data=mask_hdu['SCI'].data
    mask=np.zeros((2048,2048))
    for j in range(2048):
            ra_pix,dec_pix = wcs.pixel_to_world_values(np.arange(2048),np.zeros(2048)+j) ##Changed 2048 to 2038 to not mask borders
            x_mask,y_mask = wcs_mask.all_world2pix(ra_pix,dec_pix,0)
            if np.any((x_mask > mask_data.shape[1])|(y_mask > mask_data.shape[0])):
                print(np.max(x_mask), np.max(y_mask), mask_data.shape, j)
            y_round = np.array(np.round(y_mask),dtype='int')
            y_round[(y_round >= mask_data.shape[0])] = mask_data.shape[0]-1
            x_round = np.array(np.round(x_mask),dtype='int')
            x_round[(x_round >= mask_data.shape[1])] = mask_data.shape[1]-1
            mask[j,:]=mask_data[y_round,x_round]

    return mask == 1
#################################################################################
def med_peramp(inimg, src_mask, dq_data):
    #per amp subtraction
    img = np.copy(inimg)

    dq_mask = (build_mask(dq_data,  '~DO_NOT_USE+NON_SCIENCE')==0)
    comb_mask = (src_mask)|(dq_mask)

    bkg_img = np.zeros(img.shape, dtype=float)

    img_npix = np.zeros(img.shape, dtype=int)
    #do row sub for each amplifier
    min_frac_save = 1.
    for namp in range(4):
        mean, median, stddev = sigma_clipped_stats(img[:,512*namp:512*(namp+1)], mask=comb_mask[:,512*namp:512*(namp+1)], axis=1)
        #grow into same shape 
        bkg_img[:,512*namp:512*(namp+1)] = np.array([median]*512).T
        #Count min number of unmasked pixesls
        row_fully_masked = np.all(dq_mask[:,512*namp:512*(namp+1)], axis=1)     #rows to remove from count because all dq masked
        ngood = np.sum(( comb_mask[:,512*namp:512*(namp+1)] == False), axis=1)
        img_npix[:,512*namp:512*(namp+1)] = np.array([ngood]*512).T
        min_pix = np.nanmin(ngood[(row_fully_masked == False)])
        frac_good = ngood/np.sum(dq_mask[:,512*namp:512*(namp+1)]==False, axis=1)
        min_frac = np.nanmin(frac_good[(row_fully_masked == False)])

        if np.nanmin(frac_good) < 0.4:
            n_low = np.sum(frac_good[(row_fully_masked == False)] <0.4)
            n_vlow = np.sum(frac_good[(row_fully_masked == False)] < 0.1)
            #print('Warming, low npix ', min_pix, n_low, ' below 0.4 ', n_vlow, ' below 0.1')
        if min_frac < min_frac_save:
            min_frac_save = min_frac
        if min_pix < 100:
            fix_row =  np.arange(2048)[(ngood < 100)&(~row_fully_masked)]
            for irow in fix_row:
                mean, median, stddev = sigma_clipped_stats(img[irow,:], mask=comb_mask[irow,:])
                bkg_img[irow,512*namp:512*(namp+1)]  = median
                #print('using full row ', irow, ' amp ', namp, median)
                ngood = np.sum(( comb_mask[irow,:] == False))
                #img_npix[irow,512*namp:512*(namp+1)] = ngood
            

    #only sub unmasked
    sub_mask = (build_mask(dq_data,  '~DO_NOT_USE+NON_SCIENCE')==1) #dont apply to do not use
    sub_img_row = np.copy(img)
    #sub_img_row[(sub_mask)] = sub_img_row[(sub_mask)] - bkg_img[(sub_mask)]
    sub_img_row = sub_img_row - bkg_img

    #col sub
    mean, median, stddev = sigma_clipped_stats(sub_img_row, mask=comb_mask, axis=0)
    #grow into same shape 
    bkg_vert = np.array([median]*2048)
    sub_img_rowcol = np.copy(sub_img_row)
    #sub_img_rowcol[(sub_mask)] = sub_img_rowcol[(sub_mask)] - bkg_vert[(sub_mask)]
    sub_img_rowcol = sub_img_rowcol - bkg_vert

    bkg_amp = np.zeros(img.shape, dtype=float)
    for namp in range(4):
        mean, median, stddev = sigma_clipped_stats(sub_img_rowcol[:,512*namp:512*(namp+1)], mask=comb_mask[:,512*namp:512*(namp+1)])
        #grow into same shape 
        bkg_amp[:,512*namp:512*(namp+1)] = median
    sub_img_rowcolamp = np.copy(sub_img_rowcol)
    #sub_img_rowcolamp[(sub_mask)] = sub_img_rowcolamp[(sub_mask)] - bkg_amp[(sub_mask)]
    sub_img_rowcolamp = sub_img_rowcolamp - bkg_amp

    return sub_img_rowcolamp, bkg_img, bkg_vert, bkg_amp, img_npix
#################################################################################
def single_fix(params, save_info=False):
    filename      = params['filename']
    wisp_template = params['wisp_template']
    save_dir      = params['save_dir']
    mask_hdr      = params['mask_hdr']
    mask          = params['mask']

    cal_dm = ImageModel(filename)
    original_image = np.copy(cal_dm.data)

    #apply wisp template
    cal_dm.data -= wisp_template

    #make sball mask
    sball_mask = make_snowball_mask(cal_dm.dq)      #make_snowball_mask(cal_hdu['DQ'].data)
    cal_dm.dq[(sball_mask)] = 1                     ##cal_hdu['DQ'].data[(sball_mask)] = 1

    #filter
    
    src_mask = make_src_mask_dm(cal_dm, mask, mask_hdr)
    #presubtraction
    dq_mask = (build_mask(cal_dm.dq,  '~DO_NOT_USE+NON_SCIENCE')==0)  #(build_mask(cal_hdu['DQ'].data,  '~DO_NOT_USE+NON_SCIENCE')==0)
    mean, median, stddev = sigma_clipped_stats(cal_dm.data, mask=(dq_mask)|(src_mask), sigma=3.0)
    #cal_hdu['SCI'].data -= median
    cal_dm.data -= median
    #filter data
    #cal_hdu['SCI'].data, bkg_sub, bkg_vert, bkg_amp, img_npix = med_peramp(cal_hdu['SCI'].data, src_mask, cal_hdu['DQ'].data)
    cal_dm.data, bkg_sub, bkg_vert, bkg_amp, img_npix = med_peramp(cal_dm.data, src_mask, cal_dm.dq)

    #set median to zero
    #set median sky to zero
    #dq_mask = (build_mask(cal_dm.dq,  '~DO_NOT_USE+NON_SCIENCE')==0)  #(build_mask(cal_hdu['DQ'].data,  '~DO_NOT_USE+NON_SCIENCE')==0)
    mean, median, stddev = sigma_clipped_stats(cal_dm.data, mask=(dq_mask)|(src_mask), sigma=3.0)
    #cal_hdu['SCI'].data -= median
    cal_dm.data -= median
    #update meta data
    cal_dm.meta.background.level = 0.0
    cal_dm.meta.background.subtracted = True

    #write fixed cal file
    filtname = save_dir + (filename.split('/')[-1])
    #cal_hdu.writeto(filtname, overwrite=True)
    cal_dm.data[(original_image == 0.)] = 0.
    cal_dm.write(filtname)

    #save debug info
    if save_info:
        hdu1 = fits.PrimaryHDU(bkg_sub)
        hdu2 = fits.ImageHDU(bkg_vert)
        hdu3 = fits.ImageHDU(bkg_amp)
        hdu4 = fits.ImageHDU(img_npix)

        hdul = fits.HDUList([hdu1, hdu2, hdu3, hdu4])
        hdul.writeto(filtname.replace('crf.fits', 'fix_info.fits'), overwrite=True)
#################################################################################
def fix_cals_para(FILTER, mask_hdu, wisp_dir, save_dir, data_dir, n_procs=40):

    if FILTER in ['F277W', 'F356W', 'F410M', 'F444W']: 
        camlist = ['long']
    else:
        camlist = [1,2,3,4]

    for module in ['a','b']:
        for camera in camlist:
            #load wisp template
            if camera != 'long':
                wisp_template_raw = fits.getdata(wisp_dir + 'cal2_srcmask_median_stack_nrc%s%s.fits'%(module,camera))
                biwt = biweight_location(wisp_template_raw, ignore_nan=True)
                wisp_template = wisp_template_raw - biwt
            else:
                wisp_template = np.zeros([2048,2048], dtype=float)

            fits_files_cal = glob(data_dir+'jw*nrc%s%s_crf.fits'%(module,camera))

            with Pool(n_procs) as pool:
                dictlist = [dict(filename=fname, wisp_template=wisp_template, mask=mask_hdu['SCI'].data, mask_hdr=mask_hdu['SCI'].header, save_dir=save_dir) for fname in fits_files_cal]
                pool.map(single_fix, dictlist)
#################################################################################
def restack(params):
    #load wcs 
    FILTER = params['FILTER']
    basedir = params['basedir']
    save_dir = basedir+FILTER+'/pipe4_filt/'

    #load wcs 
    wcs_pickle = basedir + 'output_grid_wcs.npz'
    save_wcs = np.load(wcs_pickle)

    #Step 3:
    fits_files_corr = glob(save_dir + 'jw*crf.fits')

    pipe = Image3Pipeline()
    pipe.output_dir=save_dir
    pipe.tweakreg.skip=True
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

    os.system('mv %sstep_i2d.fits %sstack_%s_pipe4_v1_%s.fits'%(save_dir, save_dir, FILTER, datetime.today().strftime('%Y%m%d')))
#################################################################################   
#def dep():
    #restack('F356W')

    #test it
    # fname = '/scratch/mruari/EIGER/imaging/J0100+2802/F200W/pipe2_mywcs/jw01243001001_04101_00002_nrcb1_crf.fits'
    # save_dir = '/scratch/mruari/EIGER/imaging/J0100+2802/junk/outtest/'
    # mask_hdu = fits.open('/scratch/mruari/EIGER/imaging/J0100+2802/F200W/mycals/sourcemask_pipe3_v1.fits')
    # wisp_template = np.zeros([2048,2048], dtype=float)
    # params = dict(filename=fname, wisp_template=wisp_template, mask_hdu=mask_hdu, save_dir=save_dir) 
    # single_fix(params)

    # #check missing files
    # FILTERS = ['F115W', 'F200W', 'F356W']
    # for FILTER in FILTERS:

    #     pipe4_dir = basedir+FILTER+'/pipe4_filt/'
    #     pipe2_dir = basedir+FILTER+'/pipe2_mywcs/'
    #     wisp_dir = basedir+FILTER+'/mycals/'

    #     fits_pipe2 = glob(pipe2_dir + 'jw*crf.fits')
    #     for fname in fits_pipe2:
    #         if not os.path.isfile(pipe4_dir + fname.split('/')[-1]):
    #             print('missing file ', FILTER, fname.split('/')[-1])
                 
    #             #open mask
    #             mask_name =  wisp_dir + 'sourcemask_pipe3_v1.fits'
    #             mask_hdu = fits.open(mask_name)

    #             if FILTER == 'F356W':
    #                 wisp_template = np.zeros([2048,2048], dtype=float)
    #             else:
    #                 module = fname[fname.find('_nrc')+4]
    #                 camera = fname[fname.find('_nrc')+5]
    #                 #open wisp
    #                 wisp_template_raw = fits.getdata(wisp_dir + 'cal2_srcmask_median_stack_nrc%s%s.fits'%(module,camera))
    #                 biwt = biweight_location(wisp_template_raw, ignore_nan=True)
    #                 wisp_template = wisp_template_raw - biwt

    #             params = dict(filename=fname, wisp_template=wisp_template, mask=mask_hdu['SCI'].data, mask_hdr=mask_hdu['SCI'].header, save_dir=pipe4_dir) 
    #             single_fix(params)
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

        save_dir = basedir+FILTER+'/pipe4_filt/'
        data_dir = basedir+FILTER+'/pipe2_mywcs/'
        wisp_dir = basedir+FILTER+'/mycals/'

        if not os.path.isdir(save_dir): os.mkdir(save_dir)
        os.chdir(save_dir)

        mask_name =  wisp_dir + 'sourcemask_cal2_v1.fits'

        #open mask
        mask_hdu = fits.open(mask_name)

        #fix cal files
        fix_cals_para(FILTER, mask_hdu, wisp_dir, save_dir, data_dir)

    #check files complete
    missing = False
    for FILTER in FILTERS:

        pipe4_dir = basedir+FILTER+'/pipe4_filt/'
        pipe2_dir = basedir+FILTER+'/pipe2_mywcs/'
        wisp_dir = basedir+FILTER+'/mycals/'

        fits_pipe2 = glob(pipe2_dir + 'jw*crf.fits')
        for fname in fits_pipe2:
            if not os.path.isfile(pipe4_dir + fname.split('/')[-1]):
                print('missing file ', FILTER, fname.split('/')[-1])
                missing = True

    if missing:
        sys.exit('missing filtered files')

    #now stack
    # with Pool(3) as pool:
    #     dictlist = [dict(FILTER=FILT, basedir=basedir) for FILT in FILTERS]
    #     pool.map(restack, dictlist)
    restack(dict(FILTER='F356W', basedir=basedir))
    restack(dict(FILTER='F200W', basedir=basedir))
    restack(dict(FILTER='F115W', basedir=basedir))

##########################################################  
if __name__ == "__main__":
    main()
