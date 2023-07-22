#from astropy.visualization import simple_norm
#import matplotlib.pyplot as plt
import numpy as np
from photutils.detection import find_peaks
from astropy.nddata import NDData
from astropy.io import fits
from astropy.table import Table
from photutils.psf import extract_stars
from photutils.centroids import centroid_2dg, centroid_com
from photutils.psf import EPSFBuilder
from astropy.stats import sigma_clipped_stats
from photutils.psf import CosineBellWindow, create_matching_kernel
import os
#from astropy.visualization import LogStretch
#from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.convolution import Gaussian2DKernel, convolve
from photutils.psf.matching import resize_psf
from photutils import CircularAperture, aperture_photometry, CircularAnnulus
from photutils.psf import (HanningWindow, TukeyWindow, CosineBellWindow,
                           SplitCosineBellWindow, TopHatWindow)
#from matplotlib.gridspec import GridSpec
#from matplotlib import colors
import yaml
import argparse
import webbpsf
import sys
from glob import glob
###################################################
def get_model_JWST(FILTER, oversamp_fact=1, SNR=None, fov_arcsec=2):
    #get F356W PSF with extra fact 2 in oversamp
    instrument = webbpsf.NIRCam()
    instrument.filter=FILTER
    instrument.options['parity'] = 'odd'
    instrument.pixelscale = 0.03

    psf_model = instrument.calc_psf(oversample=oversamp_fact, fov_arcsec=fov_arcsec)

    return psf_model['DET_SAMP'].data
###################################################
def get_model_HST(basedir, FILTER):
    PSF_fits = fits.open(basedir+'HST/PSF/'+FILTER+'/'+FILTER+'_tinytim00.fits')

    PSF_model = PSF_fits[0].data                                            #0.01 arcsec
    diff_kernel = np.array([PSF_fits[0].header[-3].split(),                     #0.05 arcsec
                            PSF_fits[0].header[-2].split(),
                            PSF_fits[0].header[-1].split()]).astype(float)


    #resize to match mosaic sampling
    PSF_model_match = resize_psf(PSF_model, 0.0101, 0.03)
    diff_kernel_match = resize_psf(diff_kernel, 0.05, 0.03)

    #convolve PSF
    PSF_model_diff = convolve(PSF_model_match, diff_kernel_match, normalize_kernel=True, preserve_nan=True)
    return PSF_model_diff
###################################################
def prep_file_JWST(source_name, save_name, FILTER, ref_PSF, alpha=0.2, beta=0.6):
    img = fits.open(source_name)
    data_copy = np.copy(img['SCI'].data)

    PSF = get_model_JWST(FILTER)

    #change zeros to nans
    data_copy[(img['WHT'].data == 0)] = np.nan

    #get kernel
    window = SplitCosineBellWindow(alpha=alpha, beta=beta)
    PSF_kernel = create_matching_kernel(PSF, ref_PSF, window=window)

    #convolve
    data_conv = convolve(data_copy, PSF_kernel, normalize_kernel=True, preserve_nan=True)

    #save
    img['SCI'].data = data_conv
    img.writeto(save_name)
###################################################
def prep_file_HST(source_name, save_name, FILTER, ref_PSF, basedir,alpha=0.2, beta=0.6):
    img = fits.open(source_name)
    data_copy = np.copy(img[0].data)

    PSF_full = get_model_HST(basedir, FILTER)
    #crop to match ref_PSF
    crop_size = int((PSF_full.shape[0] - ref_PSF.shape[0])/2)
    PSF = PSF_full[crop_size:-crop_size,crop_size:-crop_size]

    #change zeros to nans
    data_copy[(img[0].data == 0)] = np.nan

    #get kernel
    window = SplitCosineBellWindow(alpha=alpha, beta=beta)
    PSF_kernel = create_matching_kernel(PSF, ref_PSF, window=window)

    #convolve
    data_conv = convolve(data_copy, PSF_kernel, normalize_kernel=True, preserve_nan=True)

    #save
    img[0].data = data_conv
    img.writeto(save_name)
##########################################################  
def write_paramfile(basedir, ver, det_band, nrc_stacks, hst_stacks):
    data =  {            'basedir': basedir,
                         'ver': ver,
                         'det_band':det_band,
                         'nrc_stacks': nrc_stacks,
                         'hst_stacks': hst_stacks}

    with open('photometry_params_v%s.yml'%(ver), 'w') as yaml_file:
        yaml.dump(data, yaml_file, default_flow_style=False)

    # with open('reduction_params.yml') as f:
    # # use safe_load instead load
    #     dataMap = yaml.safe_load(f)
########################################################## 
def main():
    p = argparse.ArgumentParser()
    p.add_argument("-v", "--ver", help='version number of catalog 1 ')
    p.add_argument("-d", "--dir", help='base directory e.g. "/scratch/mruari/EIGER/imaging/J0100+2802/"')


    args = p.parse_args() 


    basedir = args.dir
    photdir = basedir + 'photometry/'
    workdir = photdir+'PSF_matching/'
    if not os.path.isdir(workdir): os.mkdir(workdir)
    stackdir  = basedir+'junk/current_best/'
    maskdir  = stackdir+'src_mask/'

    #find stacks
    det_band = 'F356W'
    nrc_stacks = []

    for FILTER in ['F356W', 'F115W', 'F200W']:

        filenames = glob(stackdir+'stack_%s_pipe4_*.fits'%(FILTER))
        if len(filenames) != 1: 
            sys.exit('error, more than one current_best %s'%(FILTER))
        fname = filenames[0].split('/')[-1]

        if FILTER == det_band:   
            detection = True
            conv_name = stackdir + fname
        else:                   
            detection = False
            conv_name = workdir + fname.replace('.fits','.conv.fits')


        masks = glob(maskdir+'%s_sourcemask*.fits'%(FILTER))
        if len(masks) != 1: 
            sys.exit('error, more than one sourcemask %s'%(FILTER))
        mname = masks[0]

        nrc_stacks.append(dict(FILTER=FILTER,  mosaic=stackdir+fname, src_mask=mname, conv_name=conv_name, detection=detection,  inst='nrc', zpt=28.03))

    #load HST data
    #TBD
    hst_stacks = []

    for FILTER in ['F775W']:


        fname = '%s_stack_gaiaref_nrctweak_final_skylocalmin2_drc_sci.fits'%(FILTER)
        #get zpt
        hdr = fits.getheader(stackdir + fname)
        zpt = hdr['ZPT_AB']

        detection = False
        conv_name = workdir + fname.replace('.fits','.conv.fits')


        masks = glob(maskdir+'%s_sourcemask*.fits'%('F115W'))
        if len(masks) != 1: 
            sys.exit('error, more than one sourcemask %s'%('F115W'))
        mname = masks[0]

        hst_stacks.append(dict(FILTER=FILTER,  mosaic=stackdir+fname, src_mask=mname, conv_name=conv_name, detection=detection,  inst='acs', zpt=zpt))


    print( args.ver, 'photometry_params_v%s.yml'%( args.ver ) )
    if not os.path.isfile('photometry_params_v%s.yml'%( args.ver )):
        write_paramfile(args.dir, args.ver, det_band, nrc_stacks, hst_stacks)

    # #do convolutions
    psf_F356W = get_model_JWST('F356W')

    # FILTERS = ['F115W', 'F200W']
    for stackdata in nrc_stacks[1:]:
        source_name = stackdata['mosaic']
        save_name = stackdata['conv_name']
        prep_file_JWST(source_name, save_name, stackdata['FILTER'], psf_F356W)

    #FILTERS = ['F606W', 'F775W','F850LP'] 
    for stackdata in hst_stacks:
        #fname = FILTER+'_stack_F356Wref_final_drc_sci.fits'
        source_name = stackdata['mosaic']
        save_name = stackdata['conv_name']
        prep_file_HST(source_name, save_name, stackdata['FILTER'], psf_F356W, basedir)
##########################################################  
if __name__ == "__main__":
    main()