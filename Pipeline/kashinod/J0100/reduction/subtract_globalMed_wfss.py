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
from scipy.optimize import curve_fit
import argparse


## This script is for processing WFSS rat.fits or cal.fits files
## This can be used only for GRISMR.
usage = 'Usage: python subtract_globalMed_wfss.py input_rate/cal.fits out_dir --mask mask_name(optional)'
parser = argparse.ArgumentParser()


try:
    parser.add_argument("inp_fil", type=str)
    parser.add_argument("out_dir", type=str)
    parser.add_argument("--mask", type=str)
    args=parser.parse_args()
except Exception as e:
    print(usage)
    print(e)
    
inp_fil = args.inp_fil
out_dir = args.out_dir

## Output directory
if os.path.isdir(out_dir) == False:
    print('Make directory: ', out_dir, flush=True)
    os.mkdir(out_dir)
else:
    print('Output directory already exits: ', out_dir, flush=True)


## Mask image
if args.mask!=None:
    masking=True
    mask_fil = args.mask
    if os.path.isfile(mask_fil) == False:
        raise ValueError('No mask file exists: '+ mask_fil)
else:
    masking=False
    print('No mask is applied.')

## Read input fits file name
fits_name = inp_fil.split("/")[-1]
out_fil = os.path.join(out_dir, fits_name)

## Output file name
## File names for intermediate products (ip)
## Global X/Y medians

globalMed_fil = out_fil.replace(".fits", "_globalMed.fits")

## Original - global X/Y medians
globalMedSubt_fil = out_fil.replace(".fits", "_globalMedSubt.fits")

## PDF figure
## plt_fil = out_fil.replace(".fits",".pdf")

print('Input fits file : ', inp_fil, flush=True)
print('Output global_median fits file: ', globalMed_fil, flush=True)
print('Output global-median-subtracted fits file: ', globalMedSubt_fil, flush=True)

# Open fits file
# EXT 1 SCI
# EXT 2 ERR
# EXT 3 DQ
# EXT 4 VAR_POISSON
# EXT 5 VAR_RNOISE
# EXT 6 ASDF

hdul = fits.open(inp_fil)
imheader = hdul[1].header
sci_img = np.copy(hdul[1].data)
err_img = np.copy(hdul[2].data)
dq_img  = np.copy(hdul[3].data)

## Image shape
n_y, n_x = sci_img.shape[0], sci_img.shape[1]
print('n_x, n_y=', n_x, n_y,  '(', fits_name, ')', flush=True)

## Copy arrays
sci_img_original = np.copy(sci_img)
sci_img_bg = np.copy(sci_img)  ## The global (column) median will be taken from this sci_img_bg

## Read mask image
if masking:
    print('Read mask_fil: ', mask_fil)
    mask = fits.getdata(mask_fil,0)
    idx_maskedout = np.where(mask==1)  ## Pixels of bright emission lines
    idx_unmasked = np.where(mask==0)   
    print('# idx_maskedout: ', idx_maskedout[0].size)
    print('# idx_unmasked : ', idx_unmasked[0].size)
    sci_img_bg[idx_maskedout]=np.nan
else:
    idx_maskedout= np.where(np.isnan(sci_img_bg))
    idx_unmasked = np.where(np.isfinite(sci_img_bg))

## Bad pixels and any pixels with negative error
idx_badpixels = np.where((err_img <= 0)|(np.mod(dq_img, 2)==1))
print('Bad pixels: ', idx_badpixels[0].size, '(', fits_name, ')', flush=True)

## Replace bad pixels with nan
sci_img_bg[idx_badpixels]=np.nan

## Get percentiles to remove outlier pixels
err_ptiles = np.nanpercentile(err_img[idx_unmasked], [0.15, 1., 2.5, 16., 50., 84., 97.5, 99., 99.85, 99.99])
sci_ptiles = np.nanpercentile(sci_img[idx_unmasked], [0.15, 1., 2.5, 16., 50., 84., 97.5, 99., 99.85, 99.99])
print('ERR Ptiles: ', err_ptiles)
print('SCI Ptiles: ', sci_ptiles)

## Outliers if outside the 1-99th percentile
idx_outliers = np.where((sci_img<sci_ptiles[1])|(sci_img>sci_ptiles[7])|
                        (err_img<err_ptiles[1])|(err_img>err_ptiles[7]))

## Replace outliers with nan
sci_img_bg[idx_outliers]=np.nan

## Empty array for the column median
sci_img_medians_at_x = np.zeros((n_y, n_x))  # Medians at each "X", thus "vertically" striped

print(datetime.now(), ': Start - getting global median subtraction', flush=True)
for i_x in range(n_x):
    sci_img_medians_at_x[:,i_x] = np.nanmedian(sci_img_bg[:,i_x])
print(datetime.now(), ': End - getting global median subtraction', flush=True)

globalMed_img =  sci_img_medians_at_x
globalMedSubt_img = sci_img_original - sci_img_medians_at_x

### Save global (column) median subtracted image
### with all the other extensions
hdul[1].data = globalMedSubt_img
hdul.writeto(globalMedSubt_fil, overwrite=True)
print(datetime.now(), ': Saved: ', globalMedSubt_fil)

### Save global (column) median image (only; without other extensions)
hdul[1].data = globalMed_img
hdul[1].writeto(globalMed_fil, overwrite=True)
print(datetime.now(), ': Saved: ', globalMed_fil)


