from astropy.io import fits
import numpy as np
import sys, os, copy
from datetime import datetime
from astropy.convolution import interpolate_replace_nans
from astropy.convolution import convolve
from astropy.convolution import convolve_fft
from datetime import datetime
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1
from scipy.optimize import curve_fit
import argparse


### Function of median filtering for 2D array
def median_filter2d(array, kx, ky, kx_gap=0, preserve_nan=False):

    ### kx, ky: the kernel size in x (column) and y (row) direction
    ### For NIRCam GRISMR, ky should be 1.
    ### kx_gap is the size of the "hole" at the center x.
    ### preserve_nan; if True, after performing median filtering, pixels that were originally NaN again become NaN.
    
    ### Convert the input array into a numpy array
    array = np.array(array)

    ### Get shape
    n_y0, n_x0 = array.shape

    ### Add 1 pixel at the edges the input array, and fill them with NaN.
    ### NaN will be ignored
    n_y = n_y0 + 2
    n_x = n_x0 + 2
    array_ext = np.zeros((n_y,n_x))+np.nan
    array_ext[1:-1, 1:-1] = array

    ### Indices arrays
    idx_x = np.fromfunction(lambda  i, j: i+j, (n_x, kx), dtype=np.int64) - kx // 2
    idx_y = np.fromfunction(lambda  i, j: i+j, (n_y, ky), dtype=np.int64) - ky // 2

    ### Exclude the center gap
    if kx_gap > 0:
        #print('Gap size: ', kx_gap, flush=True)
        cols_remain=np.append(np.arange((kx - kx_gap)/2., dtype=np.int64),
                              np.flip(kx-1-np.arange((kx - kx_gap)/2., dtype=np.int64)))
        idx_x = idx_x[:, cols_remain]

    idx_x[idx_x < 0]=0
    idx_x[idx_x > n_x-1]=n_x-1
    idx_y[idx_y < 0]=0
    idx_y[idx_y > n_y-1]=n_y-1

    out_img = np.zeros(array_ext.shape)
    for iy in np.arange(n_y):
        med_arr = np.nanmedian(array_ext[idx_y[iy,0]:idx_y[iy,-1]+1,idx_x], axis=[0,2])
        out_img[iy,:] = med_arr
        
    if preserve_nan==True:
        out_img[np.isnan(array_ext)]=np.nan

    return out_img[1:-1,1:-1]  


## This script is for processing WFSS rat.fits or cal.fits files
## This can be used only for GRISMR.
## This should be run after subtract_globalMed_wfss.py

## usage = 'Usage: python subtract_continua_wfss.py input.fits(rate/cal) out_dir kernel_x kernel_y kernel_x_gap --mask mask_name(optional) --short_kernel kx_s, ky_s, kx_gap_s'

### Kernel parameters and those for identifying large-gradient pixels must be adjusted carefully.
### Preliminary parameters are
### kx = 51
### ky = 1
### kx_gap = 9
### Short kernel
### kx_s = 21
### ky_s = 1
### kx_gap_s = 5

### The parameters for gradient map.
### These are subject to the kernel sizes.
### 2022.08.24 Adjusted by D.Kashino using the J0100+2802 data.
gradientmap_params = {'npix_meas':50, 
                      'kx_meas':101,            ## Median kernel for measuring gradients
                      'ky_meas':1,
                      'kx_gap_meas':9,
                      'mean_SN_thresh': 3.0,     ## S/N threshold
                      'gradient_thresh': 0.3,    ## gradient threshold (per dx = 50 pixel)
                      'kernel_shape1': (1, 21),  ## For cleaning (exclude outlier pixels in the gradient map)
                      'kernel_shape2': (3, 101), ## For expanding the flagged region
                      'kernel_shape3': (3, 101)} ## For smoothing


### Get argument
parser = argparse.ArgumentParser()

try:
    parser.add_argument("inp_fil", type=str)
    parser.add_argument("out_dir", type=str)
    parser.add_argument("kx", type=int)
    parser.add_argument("ky", type=int)
    parser.add_argument("kx_gap", type=int)
    parser.add_argument("--mask", type=str)
    parser.add_argument("--short_kernel", type=int, nargs="*")
    args=parser.parse_args()
except Exception as e:
    print(e)


## Input file and output directory
inp_fil = args.inp_fil
out_dir = args.out_dir

## Output directory
if os.path.isdir(out_dir) == False:
    print('Make directory: ', out_dir, flush=True)
    os.mkdir(out_dir)
else:
    print('Output directory already exits: ', out_dir, flush=True)


## Kernel
kx = args.kx
ky = args.ky
kx_gap = args.kx_gap

## kx ky and kx_gap must be a odd number
if kx%2 == 0 or ky%2 == 0 or kx_gap%2 == 0:
    raise ValueError('kx, ky, and kx_gap all must be an odd number (kx_gap must be =>1).')
if kx <= kx_gap:
    raise ValueError('kx_gap must be smaller than kx.')

## Note for recording the parameters used for this process in the header of the output fits files
CONTSUBT_NOTE='Kernel kx {}, ky {}, kx_gap {}'.format(kx,ky,kx_gap)

## Shorter kernel for high-gradient region
if args.short_kernel:

    ## Short kernel sizes
    kx_s = args.short_kernel[0]
    ky_s = args.short_kernel[1]
    kx_gap_s = args.short_kernel[2]

    ## kx ky and kx_gap must be a odd number
    if kx_s%2 == 0 or ky_s%2 == 0 or kx_gap_s%2 == 0:
        raise ValueError('kx_s, ky_s, and kx_gap_s all must be an odd number (kx_gap_s must be =>1).')
    if kx_s <= kx_gap_s:
        raise ValueError('kx_gap_s must be smaller than kx_s.')

    CONTSUBT_NOTE=CONTSUBT_NOTE+'; Short kernel kx_s {}, ky_s {}, kx_gap_s {}'.format(kx_s,ky_s,kx_gap_s)
    CONTSUBT_NOTE=CONTSUBT_NOTE+' with gradient map params - '\
        'npix_meas {}, '\
        'kx_meas {}, '\
        'ky_meas {}, '\
        'kx_gap_meas {}, '\
        'mean_SN_thresh {}, '\
        'gradient_thresh {}, '\
        'kernel_shape1 {}, '\
        'kernel_shape2 {}, '\
        'kernel_shape3 {}'.format(gradientmap_params['npix_meas'],
                                  gradientmap_params['kx_meas'],
                                  gradientmap_params['ky_meas'],
                                  gradientmap_params['kx_gap_meas'],
                                  gradientmap_params['mean_SN_thresh'],
                                  gradientmap_params['gradient_thresh'],
                                  gradientmap_params['kernel_shape1'],
                                  gradientmap_params['kernel_shape2'],
                                  gradientmap_params['kernel_shape3'])
                                                                                                                                                                               
    print('Apply short kernel for pixels which have large gradient: ', kx_s, ky_s, kx_gap_s, flush=True)
    print('Short kernel (kx, ky, kx_gap): ', kx_s, ky_s, kx_gap_s, flush=True)



## Mask image
if args.mask!=None:
    masking=True
    mask_fil = args.mask
    if os.path.isfile(mask_fil) == False:
        raise ValueError('No mask file exists: '+ mask_fil)
else:
    masking=False
    print('No mask is applied.')


    
## Input fits file name
fits_name = inp_fil.split("/")[-1]



## Output file name
out_fil = os.path.join(out_dir, fits_name)
emline_fil = out_fil.replace(".fits", "_emline.fits")

## File name for "continua image"
continua_fil = out_fil.replace(".fits", "_continua.fits")

print('Input fits file                       : ', inp_fil, flush=True)
print('Output emline   fits file             : ', emline_fil, flush=True)
print('Output continua fits file             : ', continua_fil, flush=True)


## File name for "continua image with long kernel"
if args.short_kernel:
    continua_wht_fil = out_fil.replace(".fits", "_continua_wht.fits")
    continua_lk_fil = out_fil.replace(".fits", "_continua_lk.fits")
    continua_sk_fil = out_fil.replace(".fits", "_continua_sk.fits")
    print('Output continua weight fits file      : ', continua_wht_fil, flush=True)
    print('Output continua long-kernel fits file : ', continua_lk_fil, flush=True)
    print('Output continua short-kernel fits file: ', continua_sk_fil, flush=True)


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

## SCI background (bkg) image
## The continua image will be created from this bkg_img
bkg_img = copy.deepcopy(sci_img)

## Image shape
n_y, n_x = sci_img.shape[0], sci_img.shape[1]
print('n_x, n_y=', n_x, n_y,  '(', fits_name, ')', flush=True)


## Read mask image
if masking:
    print('Read mask_fil: ', mask_fil)
    mask = fits.getdata(mask_fil,0)
    idx_maskedout = np.where(mask==1)  ## Pixels of bright emission lines
    idx_unmasked = np.where(mask==0)   
    print('# idx_maskedout: ', idx_maskedout[0].size)
    print('# idx_unmasked : ', idx_unmasked[0].size)
    bkg_img[idx_maskedout]=np.nan

## Bad pixels and any pixels with negative error
idx_badpixels = np.where((err_img <= 0)|(np.mod(dq_img, 2)==1))
print('Bad pixels: ', idx_badpixels[0].size, '(', fits_name, ')', flush=True)

## Replace bad pixels with nan
sci_img[idx_badpixels]=np.nan
bkg_img[idx_badpixels]=np.nan
err_img[idx_badpixels]=np.nan

print(datetime.now(), ': Start continuum subtraction')

### Continuum map using median_filter2d
continua_img = median_filter2d(bkg_img, kx, ky, kx_gap, preserve_nan=False)

### Secondary short kernel is applied
if args.short_kernel:
    
    continua_longKernel_img = copy.deepcopy(continua_img)
    continua_shortKernel_img = median_filter2d(bkg_img, kx_s, ky_s, kx_gap_s, preserve_nan=False)

    ## Smoothed ERR image with the long kernel
    bkg_img_smoothed = median_filter2d(bkg_img,
                                       gradientmap_params['kx_meas'],
                                       gradientmap_params['ky_meas'],
                                       gradientmap_params['kx_gap_meas'],
                                       preserve_nan=False)
    err_img_smoothed = median_filter2d(err_img,
                                       gradientmap_params['kx_meas'],
                                       gradientmap_params['ky_meas'],
                                       1, preserve_nan=False)

    ## Gradient
    gradient_map = np.zeros(sci_img.shape)
    mean_map = np.zeros(sci_img.shape)

    for ix in range(sci_img.shape[1]):
        x0 = np.max([   0, ix - gradientmap_params['npix_meas']//2])
        x1 = np.min([2047, ix + gradientmap_params['npix_meas']//2])
        gradient_map[:,ix] = (bkg_img_smoothed[:,x1] - bkg_img_smoothed[:,x0]) / gradientmap_params['npix_meas']
        mean_map[:,ix]     = (bkg_img_smoothed[:,x1] + bkg_img_smoothed[:,x0]) / 2.

    ## Offset correction
    mean_map = mean_map - np.nanmin(mean_map)

    mean_map_thresholded=copy.deepcopy(mean_map)
    mean_map_thresholded[mean_map/err_img_smoothed < gradientmap_params['mean_SN_thresh']] = np.infty
    relative_gradient_map = gradient_map/mean_map_thresholded * 50 ## Don't change this 50.

    weightmap = np.zeros(relative_gradient_map.shape)
    weightmap[np.abs(relative_gradient_map) > gradientmap_params['gradient_thresh']]=1.0

    ## Removing outlier pixels (shot-noise-like pixels) 
    kernel1 = np.ones((gradientmap_params['kernel_shape1']))
    weightmap = convolve(weightmap, kernel1)
    weightmap[weightmap >= 0.5]=1.0
    weightmap[weightmap < 0.5]=0.0
    
    ## Expanding the flagged region
    kernel2 = np.ones((gradientmap_params['kernel_shape2']))
    weightmap = convolve(weightmap, kernel2)
    weightmap[weightmap >= 1e-10]=1.0
    weightmap[weightmap < 1e-10]=0.0

    ## Smoothing the boundaries
    kernel3 = np.ones((gradientmap_params['kernel_shape3']))
    weightmap = convolve(weightmap, kernel3)

    ## New continua image
    continua_img = (1 - weightmap)*continua_longKernel_img + weightmap * continua_shortKernel_img
    isnan=np.isnan(continua_shortKernel_img)
    continua_img[isnan]=continua_longKernel_img[isnan]
    
# Emission-line image (continuum extracted)
emline_img = sci_img - continua_img

# Save files
hdul[1].data = emline_img
hdul[1].header['MYNOTE']=CONTSUBT_NOTE

hdul.writeto(emline_fil, overwrite=True)
print(datetime.now(), ': Saved: ', emline_fil)

hdul[1].data = continua_img
hdul.writeto(continua_fil, overwrite=True)
print(datetime.now(), ': Saved: ', continua_fil)

if args.short_kernel:
    hdul[1].data = weightmap
    hdul.writeto(continua_wht_fil, overwrite=True)
    print(datetime.now(), ': Saved: ', continua_wht_fil)

    hdul[1].data = continua_longKernel_img
    hdul[1].writeto(continua_lk_fil, overwrite=True)
    print(datetime.now(), ': Saved: ', continua_lk_fil)

    hdul[1].data = continua_shortKernel_img
    hdul[1].writeto(continua_sk_fil, overwrite=True)
    print(datetime.now(), ': Saved: ', continua_sk_fil)

hdul.close()


