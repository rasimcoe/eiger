import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from photutils.centroids import centroid_2dg, centroid_com
from astropy.stats import sigma_clipped_stats
import os
from astropy.convolution import Gaussian2DKernel, convolve
from photutils import CircularAperture, aperture_photometry, CircularAnnulus
from scipy import interpolate
from scipy.stats import binned_statistic
import astropy.wcs as wcs
from astropy.table import QTable, Table
import astropy.units as u
from regions import Regions
from scipy.optimize import curve_fit
import yaml
import sys
import argparse
from glob import glob 
###################################################
def noise_test(imgdata, mask, radius=3.0, n_samp=1000):
    pix_grid = np.ones(imgdata.shape)

    pos_rand = np.random.rand(n_samp, 2) * np.flip(imgdata.shape)
    apers = CircularAperture(pos_rand, r=radius)
    phot_table = aperture_photometry(imgdata, apers, mask=mask)
    pix_table = aperture_photometry(pix_grid, apers, mask=mask)
    flux_samps = phot_table['aperture_sum'].data
    pixs_samps = pix_table['aperture_sum'].data

    #correct for masked area
    area_geo = np.pi * radius**2
    flux_samps = flux_samps / pixs_samps * area_geo
    #select apertures with >0.9 unmasked 
    flux_samps_select = flux_samps[(pixs_samps > area_geo*0.9)]
    clip_flux_mean, clip_flux_med, clip_flux_std = sigma_clipped_stats(flux_samps_select)
    return flux_samps, area_geo,  clip_flux_std
########################################################
def noise_funct(npix, alpha, beta):
    return alpha*(npix**beta)
########################################################
def fit_powerlaw(npix, std):
    popt, pcov = curve_fit(noise_funct, npix, std, p0=(1.0, 0.5))
    return popt, noise_funct(npix, *popt)
########################################################## 
def main():
    paramnames = glob('photometry_params_v*.yml')
    if len(paramnames) > 1: 
        sys.exit('error, more than one photometry_params_v*.yml')
    if len(paramnames) == 0: 
        sys.exit('no photometry_params_v*.yml file not found')

    with open(paramnames[0]) as f:
    # use safe_load instead load
        dataMap = yaml.safe_load(f)


    basedir = dataMap['basedir']

    nrc_stacks = dataMap['nrc_stacks']
    hst_stacks = dataMap['hst_stacks']
    basedir    = dataMap['basedir']


    workdir = basedir + 'photometry/noise_model/'
    if not os.path.isdir(workdir): os.mkdir(workdir)
    stackdir  = basedir+'junk/current_best/'
    convdir = basedir + 'photometry/PSF_matching/'
    sexdir = basedir + 'photometry/sextractor/'


    data_tests = nrc_stacks + hst_stacks       
    #random data
    radius_bins = np.logspace(0,0.6,6)
    save_flux_std = np.zeros([len(data_tests), len(radius_bins)], dtype=float)
    save_popt    =np.zeros([len(data_tests), 2], dtype=float)

    for it, test in enumerate(data_tests):
        src_mask_test = fits.getdata(test['src_mask'])
        mosaic_test = fits.open(test['conv_name'])
        if test['inst'] == 'nrc':
            #open arays
            data_copy = np.copy(mosaic_test['SCI'].data)
            #make mask
            mask = (np.isnan(data_copy))|(src_mask_test != 0)|(mosaic_test['WHT'].data == 0.)
        else:
            data_copy = np.copy(mosaic_test[0].data)
            mask = (np.isnan(data_copy))|(src_mask_test != 0)

        #bkg sub
        bkg = fits.getdata(sexdir + 'bkg_'+test['FILTER']+'.fits')
        data_copy = data_copy - bkg

        for ir, radius in enumerate(radius_bins):
            print(it, ir)
            flux_samps, npix, samp_flux_std = noise_test(data_copy, mask, radius=radius, n_samp=100000)
            save_flux_std[it, ir] = samp_flux_std

        #do fit
        popt, fit_std = fit_powerlaw(np.pi * radius_bins**2,  save_flux_std[it,:])
        save_popt[it,:] = popt 
        print(test['FILTER'], save_popt[it,:])

    filter_list = [test['FILTER'] for  test in data_tests]
    np.savez(workdir + 'noise_model_params_v%s.npz'%(dataMap['ver']), data_tests=data_tests, filters=filter_list, popt=save_popt)

##########################################################  
if __name__ == "__main__":
    main()

