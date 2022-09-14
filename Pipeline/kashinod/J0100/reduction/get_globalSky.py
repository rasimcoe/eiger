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
import argparse

# This script is for getting the median image from a list of Image2 cal.fits
# Usage: python get_globalsky.py list_of_input_files output_name [option mask_list] 

# List of cal.fits files
#list_fil = sys.argv[1]
#outname = sys.argv[2]
parser = argparse.ArgumentParser()
try:
    parser.add_argument("list_fil", type=str)
    parser.add_argument("outname", type=str)
    parser.add_argument("--mask_brightpixels", type=bool)
    parser.add_argument("--mask_list", type=str)
    args=parser.parse_args()
except Exception as e:
    print(e)

list_fil = args.list_fil
outname = args.outname

fits_files = np.loadtxt(list_fil, dtype="str")

if args.mask_list!=None: #len(sys.argv)==5:
    masking=True
    mask_list_fil = args.mask_list #sys.argv[3]
    print(datetime.now(), '-- Mask list provided: ', mask_list_fil)
    mask_fits_files = np.loadtxt(mask_list_fil, dtype="str")
    if mask_fits_files.size != fits_files.size:
        raise ValueError('Number of mask fits files must be the same as that of input fits files.')
else:
    masking=False



# Read the first file
hdul = fits.open(fits_files[0])
sci_img = np.copy(hdul[1].data)

# Prepare empty "cube"
cube = np.zeros((sci_img.shape[0], sci_img.shape[1], fits_files.size))

print(datetime.now(), '-- Start big loop.', flush=True)
for i in range(0, fits_files.size):
    print('---------------------------------------', flush=True)
    fil = fits_files[i]
    print(datetime.now(), '-- Read ', fil, flush=True)
    sci_img = fits.getdata(fil, 1)
    err_img = fits.getdata(fil, 2)

    n_y, n_x = sci_img.shape
    print(datetime.now(), '-- n_x, n_y=', n_x, n_y, flush=True)

    sci_img_original = np.copy(sci_img)
    sci_img_bg = np.copy(sci_img)

    idx_negative_err = np.where(err_img <= 0)
    print(datetime.now(), '-- Pixels with negative error: ', idx_negative_err[0].size, flush=True)
    sci_img_bg[idx_negative_err]=np.nan

    if masking:
        print(datetime.now(), '-- Read mask ', mask_fits_files[i], flush=True)
        mask = fits.getdata(mask_fits_files[i], 0)
        idx_maskedout = np.where(mask==1)
        print(datetime.now(), '-- Masked out pixels ', idx_maskedout[0].size, flush=True)
        sci_img[idx_maskedout]=np.nan
        err_img[idx_maskedout]=np.nan
        sci_img_bg[idx_maskedout]=np.nan

    if args.mask_brightpixels==True:
        err_ptiles = np.nanpercentile(err_img[:,:], [0.15, 1., 2.5, 16., 50., 84., 97.5, 99., 99.85])
        sci_ptiles = np.nanpercentile(sci_img[:,:], [0.15, 1., 2.5, 16., 50., 84., 97.5, 99., 99.85])
        print(datetime.now(), '-- ERR Ptiles: ', err_ptiles)
        print(datetime.now(), '-- SCI Ptiles: ', sci_ptiles)

        idx_bad = np.where((sci_img<sci_ptiles[1])|(sci_img>sci_ptiles[7])|
                           (err_img<err_ptiles[1])|(err_img>err_ptiles[7]))
        print(datetime.now(), '-- idx_bad[0].size: ', idx_bad[0].size)
        sci_img_bg[idx_bad]=np.nan

    cube[:,:,i]=sci_img_bg

print(datetime.now(), '-- All frames have been read.', flush=True)
print(datetime.now(), '-- Now taking medians.', flush=True)
medimg = np.nanmedian(cube, axis=2)
medimg[np.where(np.isnan(medimg))]=0.
hdul[1].data = medimg

outfil = 'globalsky_'+outname+'.fits'
print('Save '+outfil, flush=True)
hdul.writeto(outfil, overwrite=True)
