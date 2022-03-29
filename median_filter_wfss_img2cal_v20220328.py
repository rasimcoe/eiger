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


# This script is for processing WFSS Image2 cal.fits files
# This can be used only for R-grism, without considering any possible tilt.
usage = 'Usage: python median_filter_wfss_img2cas_v*.py input_cal.fits out_dir kernel_x kernel_y kernel_x_gap --mask mask_name(optional) --globalsky globalsky_name(optional)'
parser = argparse.ArgumentParser()

try:
    parser.add_argument("inp_fil", type=str)
    parser.add_argument("out_dir", type=str)
    parser.add_argument("kx", type=int)
    parser.add_argument("ky", type=int)
    parser.add_argument("kx_gap", type=int)
    parser.add_argument("--mask", type=str)
    parser.add_argument("--globalsky", type=str)
    args=parser.parse_args()
except Exception as e:
    print(e)
    

inp_fil = args.inp_fil
out_dir = args.out_dir

# kernal for median-filtering
kx = args.kx
ky = args.ky
kx_gap = args.kx_gap

# Globalsky
if args.globalsky!=None:
    subtract_globalsky=True
    globalsky_fil = args.globalsky
    if os.path.isfile(globalsky_fil) == False:
        raise ValueError('No globalsky file exists: '+ globalsky_fil)
else:
    subtract_globalsky=False
    print('No globalsky is subtracted.')


# Mask image
if args.mask!=None:
    masking=True
    mask_fil = args.mask
    if os.path.isfile(mask_fil) == False:
        raise ValueError('No mask file exists: '+ mask_fil)
else:
    masking=False
    print('No mask is applied.')

# kx ky and kx_gap must be a odd number
if kx%2 == 0 or ky%2 == 0 or kx_gap%2 == 0:
    raise ValueError('kx, ky, and kx_gap all must be an odd number.')
if kx <= kx_gap:
    raise ValueError('kx_gap must be smaller than kx.')

# Input fits file name
fits_name = inp_fil.split("/")[-1]

# Output directory
if os.path.isdir(out_dir) == False:
    print('Make directory: ', out_dir, flush=True)
    os.mkdir(out_dir)
else:
    print('Output directory already exits: ', out_dir, flush=True)


# Output file name
# The final product whre the "global x/y-medians" and "continua image" are subtracted.
# The median-subtracted file has the same name as the input.
emlines_fil = os.path.join(out_dir, fits_name)

# File name for "continua image"
continua_fil = emlines_fil.replace(".fits", "_continua.fits")

# File names for intermediate products (ip)
# Original - global X/Y medians
globalMed_fil = emlines_fil.replace(".fits", "_globalMed.fits")

# Global X/Y medians
globalMedSubt_fil = emlines_fil.replace(".fits", "_globalMedSubt.fits")

# PDF figure
plt_fil = emlines_fil.replace(".fits",".pdf")

print('Input fits file : ', inp_fil, flush=True)
print('Output global_median fits file: ', globalMed_fil, flush=True)
print('Output global-median-subtracted fits file: ', globalMedSubt_fil, flush=True)
print('Output continua fits file (original - final): ', continua_fil, flush=True)
print('Output continua-subtracted fits file (final product): ', emlines_fil, flush=True)

if inp_fil == emlines_fil:
    raise ValueError('Output file name must be different from the input.  Change the different output directory (', fits_name, ')', flush=True)

# Open fits file
# EXT 1 SCI
# EXT 2 ERR
# EXT 3 DQ
# EXT 4 VAR_POISSON
# EXT 5 VAR_RNOISE

hdul = fits.open(inp_fil)
imheader = hdul[1].header
sci_img = np.copy(hdul[1].data)
err_img = np.copy(hdul[2].data)

# Read globalsky image
if subtract_globalsky:
    print('Read globalsky: ', globalsky_fil)
    globalsky = fits.getdata(globalsky_fil,0)
    print('globalsky.shape: ', globalsky.shape)
    sci_img = sci_img[:,:] - globalsky[:,:]

n_y, n_x = sci_img.shape[0], sci_img.shape[1]
print('n_x, n_y=', n_x, n_y,  '(', fits_name, ')', flush=True)

sci_img_original = np.copy(sci_img)
sci_img_bg = np.copy(sci_img)


# Read mask image
if masking:
    mask = fits.getdata(mask_fil,0)
    idx_maskedout = np.where(mask==1)  ## Pixels of bright emission lines
    idx_unmasked = np.where(mask==0)   
    sci_img_bg[idx_maskedout]=np.nan
else:
    idx_unmasked = np.where(np.finite(sci_img_bg))
#sci_img_bg[ 0,:]=np.nan
#sci_img_bg[-1,:]=np.nan
#sci_img_bg[:, 0]=np.nan
#sci_img_bg[:,-1]=np.nan

idx_negative_err = np.where(err_img <= 0)
print('Pixels negative ERR=0: ', idx_negative_err[0].size, '(', fits_name, ')', flush=True)
sci_img_bg[idx_negative_err]=np.nan

err_ptiles = np.nanpercentile(err_img[idx_unmasked], [0.15, 1., 2.5, 16., 50., 84., 97.5, 99., 99.85, 99.99])
sci_ptiles = np.nanpercentile(sci_img[idx_unmasked], [0.15, 1., 2.5, 16., 50., 84., 97.5, 99., 99.85, 99.99])
print('ERR Ptiles: ', err_ptiles)
print('SCI Ptiles: ', sci_ptiles)

# Bad pixels if outside the 1-99th percentile
idx_bad = np.where((sci_img<sci_ptiles[1])|(sci_img>sci_ptiles[7])|
                   (err_img<err_ptiles[1])|(err_img>err_ptiles[7]))

sci_img_bg[idx_bad]=np.nan

sci_img_medians_at_x = np.zeros((n_y, n_x))  # Medians at each "X", thus "vertically" striped
sci_img_medians_at_y = np.zeros((n_y, n_x))  # Medians at each "Y", thus "horizontally" striped

print(datetime.now(), ': Start global x,y median', flush=True)

# for i_y in range(n_y):
#     sci_img_medians_at_y[i_y,:] = np.nanmedian(sci_img_bg[i_y,:])

sci_img_bg_tmp = sci_img_bg - sci_img_medians_at_y

for i_x in range(n_x):
    sci_img_medians_at_x[:,i_x] = np.nanmedian(sci_img_bg_tmp[:,i_x])

print(datetime.now(), ': End global x,y median', flush=True)

sci_img_cleaned = sci_img_original - sci_img_medians_at_x# - sci_img_medians_at_y
medians_img =  sci_img_medians_at_x# + sci_img_medians_at_y


hdu_tmp = fits.ImageHDU(data = sci_img_cleaned, header = imheader)
hdu_tmp.writeto(globalMedSubt_fil, overwrite=True)
#hdul[1].data = sci_img_cleaned
#hdul[1].writeto(globalMedSubt_fil, overwrite=True)
print(datetime.now(), ': Saved: ', globalMedSubt_fil)

hdul[1].data = medians_img
hdul[1].writeto(globalMed_fil, overwrite=True)
print(datetime.now(), ': Saved: ', globalMed_fil)


# -------------------------------------------
# Secondary "continuum" subtraction
# Kernel size
print(datetime.now(), ': Start secondary continuum subtraction')
sci_img_new_bg = np.copy(sci_img_cleaned)
sci_img_new_bg[idx_negative_err]=np.nan

idx_bad_new = np.where((sci_img<sci_ptiles[1])|(sci_img>sci_ptiles[9])|
                       (err_img<err_ptiles[1])|(err_img>err_ptiles[9]))
sci_img_new_bg[idx_bad_new] = np.nan

# Empty array for continua image
continua_img = np.zeros((n_y, n_x))

# Set the kernel
kernel=str(kx)+'x'+str(ky)
print('Kernel size: ', kernel, flush=True)

idx_x = np.fromfunction(lambda  i, j: i+j, (n_x, kx), dtype=np.int64) - kx // 2
idx_y = np.fromfunction(lambda  i, j: i+j, (n_y, ky), dtype=np.int64) - ky // 2

# Exclude the center gap
if kx_gap>0:
    print('Gap size: ', kx_gap, flush=True)
    cols_remain=np.append(np.arange((kx - kx_gap)/2., dtype=np.int64),
                           np.flip(kx-1-np.arange((kx - kx_gap)/2., dtype=np.int64)))
    idx_x = idx_x[:, cols_remain]

idx_x[idx_x < 0]=0
idx_x[idx_x > n_x-1]=n_x-1
idx_y[idx_y < 0]=0
idx_y[idx_y > n_y-1]=n_y-1

print(datetime.now(), ': Loop start (', fits_name, ')', flush=True)
for iy in np.arange(n_y):
    #print('Loop: ', iy, ' /', n_y)
    med_tmp = np.nanmedian(sci_img_new_bg[idx_y[iy,0]:idx_y[iy,-1]+1,idx_x], axis=[0,2])
    continua_img[iy,:] = med_tmp[:]
print(datetime.now(), ': Loop end (', fits_name, ')', flush=True)

# Emission-line image (continuum extracted)
sci_img_emlines = sci_img_cleaned - continua_img

# Save files
hdul[1].data = sci_img_emlines
hdul.writeto(emlines_fil, overwrite=True)
print(datetime.now(), ': Saved: ', emlines_fil)

hdul[1].data = continua_img
hdul[1].writeto(continua_fil, overwrite=True)
print(datetime.now(), ': Saved: ', continua_fil)

hdul.close()


# Plot
fig = plt.figure(figsize=(10,15))
ax = fig.subplots(3,2)
ax[0,0].set_title('Original')
vmin=sci_ptiles[1]
vmid=sci_ptiles[4]
vmax=sci_ptiles[7]

im1=ax[0,0].imshow(sci_img_original, vmin=vmin, vmax=vmax, origin="lower")

ax[0,1].set_title('Background')
kernel=np.ones((3,3))
im2=ax[0,1].imshow(interpolate_replace_nans(sci_img_bg, kernel), vmin=vmin, vmax=vmax, origin="lower")

ax[1,0].set_title('Global medians')
im3=ax[1,0].imshow(medians_img, origin="lower")


ax[1,1].set_title('Global-medians-subtracted')
im4=ax[1,1].imshow(sci_img_cleaned, vmin=vmin-vmid, vmax=vmax-vmid, origin="lower")

fig.colorbar(im1, cax=mpl_toolkits.axes_grid1.make_axes_locatable(ax[0,0]).append_axes('right','5%',pad='3%'))
fig.colorbar(im2, cax=mpl_toolkits.axes_grid1.make_axes_locatable(ax[0,1]).append_axes('right','5%',pad='3%'))
fig.colorbar(im3, cax=mpl_toolkits.axes_grid1.make_axes_locatable(ax[1,0]).append_axes('right','5%',pad='3%'))
fig.colorbar(im4, cax=mpl_toolkits.axes_grid1.make_axes_locatable(ax[1,1]).append_axes('right','5%',pad='3%'))


ax[2,0].set_title('Continua')
im5=ax[2,0].imshow(continua_img, vmin=vmin-vmid, vmax=vmax-vmid, origin="lower")

ax[2,1].set_title('emission line image')
im6=ax[2,1].imshow(sci_img_emlines, vmin=vmin-vmid, vmax=vmax-vmid, origin="lower")

fig.colorbar(im5, cax=mpl_toolkits.axes_grid1.make_axes_locatable(ax[0,0]).append_axes('right','5%',pad='3%'))
fig.colorbar(im6, cax=mpl_toolkits.axes_grid1.make_axes_locatable(ax[0,1]).append_axes('right','5%',pad='3%'))

fig.savefig(plt_fil)














