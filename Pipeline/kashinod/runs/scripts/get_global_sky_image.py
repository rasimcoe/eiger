from astropy.io import fits
import numpy as np
import sys, os
from datetime import datetime
from astropy.convolution import interpolate_replace_nans
from astropy.convolution import convolve
from astropy.convolution import convolve_fft
from datetime import datetime
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1

# This script is for getting the median image from a list of Image2 cal.fits

# List of cal.fits files
list_fil = sys.argv[1]
outname = sys.argv[2]

fits_files = np.loadtxt(list_fil, dtype="str")

do_masking=False
if len(sys.argv)==4:
    do_masking=True
    mask_list_fil = sys.argv[3]
    print('Mask list provided: ', mask_list_fil)
    mask_fits_files = np.loadtxt(mask_list_fil, dtype="str")
    if mask_fits_files.size != fits_files.size:
        raise ValueError('Number of mask fits files must be the same as that of input fits files.')

hdul = fits.open(fits_files[0])
sci_img = np.copy(hdul[1].data)
cube = np.zeros((sci_img.shape[0], sci_img.shape[1], fits_files.size))

for i in range(0, fits_files.size):
    print('-----------------------------------', flush=True)
    fil = fits_files[i]
    print('Read ', fil, flush=True)
    sci_img = fits.getdata(fil, 1)
    err_img = fits.getdata(fil, 2)

    n_y, n_x = sci_img.shape[0], sci_img.shape[1]
    print('n_x, n_y=', n_x, n_y, flush=True)

    sci_img_original = np.copy(sci_img)
    sci_img_bg = np.copy(sci_img)

    idx_negative_err = np.where(err_img <= 0)
    print('Pixels negative ERR=0: ', idx_negative_err[0].size, flush=True)
    sci_img_bg[idx_negative_err]=np.nan

    if do_masking:
        mask = fits.getdata(mask_fits_files[i], 0)
        idx_maskedout = np.where(mask==1)
        sci_img[idx_maskedout]=np.nan
        err_img[idx_maskedout]=np.nan
        sci_img_bg[idx_maskedout]=np.nan

    err_ptiles = np.nanpercentile(err_img[:,:], [0.15, 1., 2.5, 16., 50., 84., 97.5, 99., 99.85])
    sci_ptiles = np.nanpercentile(sci_img[:,:], [0.15, 1., 2.5, 16., 50., 84., 97.5, 99., 99.85])
    print('ERR Ptiles: ', err_ptiles)
    print('SCI Ptiles: ', sci_ptiles)

    idx_bad = np.where((sci_img<sci_ptiles[1])|(sci_img>sci_ptiles[7])|
                       (err_img<err_ptiles[1])|(err_img>err_ptiles[7]))

    sci_img_bg[idx_bad]=np.nan

    cube[:,:,i]=sci_img_bg


medimg = np.nanmedian(cube, axis=2)
medimg[np.where(np.isnan(medimg))]=0.
hdul[1].data = medimg
hdul.writeto('median_image_'+outname+'.fits', overwrite=True)


