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

# This script is for processing Image2 cal.fits files (after flat-fielding).
# You can provide a "global sky' image and a "mask" image.
usage = 'Usage: python median_filter_wfss_img2cas_v*.py input_cal.fits out_dir --mask mask_name(optional) --globalsky globalsky_name(optional)'
parser = argparse.ArgumentParser()

try:
    parser.add_argument("inp_fil", type=str)
    parser.add_argument("out_dir", type=str)
    parser.add_argument("--mask", type=str)
    parser.add_argument("--globalsky", type=str)
    args=parser.parse_args()
except Exception as e:
    print(e)


inp_fil = args.inp_fil
out_dir = args.out_dir

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

# Input fits file
fits_name = inp_fil.split("/")[-1]

# Output directory
out_dir = sys.argv[2]
if os.path.isdir(out_dir) == False:
    print('Make directory: ', out_dir, flush=True)
    os.mkdir(out_dir)
else:
    print('Output directory already exits: ', out_dir, flush=True)



# Output file name
# The median-subtracted file has the same name as the input 
# for subsequence reduction processes
medSubt_fil = os.path.join(out_dir, fits_name)
median_fil = medSubt_fil.replace(".fits", "_median.fits")
plt_fil = medSubt_fil.replace(".fits",".pdf")

print('Input fits file:                    ', inp_fil, flush=True)
print('Output median-subtracted fits file: ', medSubt_fil, flush=True)
print('Output median fits file:            ', median_fil, flush=True)

if inp_fil == medSubt_fil:
    raise ValueError('Output file name must be different from the input.  Change the different output directory ('+fits_name+')')


# Open fits file
# EXT 1 SCI
# EXT 2 ERR
# EXT 3 DQ
# EXT 4 VAR_POISSON
# EXT 5 VAR_RNOISE

hdul = fits.open(inp_fil)
sci_img = np.copy(hdul[1].data)
err_img = np.copy(hdul[2].data)

# Subtract globalsky image
if subtract_globalsky:
    print('Read globalsky: ', globalsky_fil)
    globalsky = fits.getdata(globalsky_fil,0)
    print('globalsky.shape: ', globalsky.shape)
    sci_img = sci_img[:,:] - globalsky[:,:]

n_y, n_x = sci_img.shape[0], sci_img.shape[1]
print('n_x, n_y=', n_x, n_y,  '(', fits_name, ')', flush=True)

sci_img_original = np.copy(sci_img)
sci_img_bg = np.copy(sci_img)

# Masking
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
print('Pixels negative error: ', idx_negative_err[0].size, '(', fits_name, ')', flush=True)
sci_img_bg[idx_negative_err]=np.nan

err_ptiles = np.nanpercentile(err_img[idx_unmasked], [0.15, 1., 2.5, 16., 50., 84., 97.5, 99., 99.85])
sci_ptiles = np.nanpercentile(sci_img[idx_unmasked], [0.15, 1., 2.5, 16., 50., 84., 97.5, 99., 99.85])
print('ERR Ptiles: ', err_ptiles)
print('SCI Ptiles: ', sci_ptiles)

idx_bad = np.where((sci_img<sci_ptiles[1])|(sci_img>sci_ptiles[7])|
                   (err_img<err_ptiles[1])|(err_img>err_ptiles[7]))

sci_img_bg[idx_bad]=np.nan

sci_img_medians_at_x = np.zeros((n_y, n_x))  # Medians at each "X", thus "vertically" striped
sci_img_medians_at_y = np.zeros((n_y, n_x))  # Medians at each "Y", thus "horizontally" striped
print(datetime.now(), '- Start global x,y median', flush=True)

for i_y in range(n_y):
    sci_img_medians_at_y[i_y,:] = np.nanmedian(sci_img_bg[i_y,:])

sci_img_bg_tmp = sci_img_bg - sci_img_medians_at_y

for i_x in range(n_x):
    sci_img_medians_at_x[:,i_x] = np.nanmedian(sci_img_bg_tmp[:,i_x])

sci_img_cleaned = sci_img_original - sci_img_medians_at_y - sci_img_medians_at_x
medians_img =  sci_img_medians_at_y + sci_img_medians_at_x


fig = plt.figure(figsize=(20,20))
ax = fig.subplots(2,2)
ax[0,0].set_title('Only global sky subtracted)')
vmin=sci_ptiles[0]
vmax=sci_ptiles[-1]
im1=ax[0,0].imshow(sci_img_original, vmin=vmin, vmax=vmax)

ax[0,1].set_title('background')
kernel=np.ones((3,3))
im2=ax[0,1].imshow(interpolate_replace_nans(sci_img_bg, kernel), vmin=vmin, vmax=vmax)

ax[1,0].set_title('medians-subtracted')
im3=ax[1,0].imshow(sci_img_cleaned, vmin=vmin, vmax=vmax)

ax[1,1].set_title('medians image')
im4=ax[1,1].imshow(medians_img)


fig.colorbar(im1, cax=mpl_toolkits.axes_grid1.make_axes_locatable(ax[0,0]).append_axes('right','5%',pad='3%'))
fig.colorbar(im2, cax=mpl_toolkits.axes_grid1.make_axes_locatable(ax[0,1]).append_axes('right','5%',pad='3%'))
fig.colorbar(im3, cax=mpl_toolkits.axes_grid1.make_axes_locatable(ax[1,0]).append_axes('right','5%',pad='3%'))
fig.colorbar(im4, cax=mpl_toolkits.axes_grid1.make_axes_locatable(ax[1,1]).append_axes('right','5%',pad='3%'))

fig.savefig(plt_fil)

hdul[1].data = sci_img_cleaned
hdul.writeto(medSubt_fil, overwrite=True)
print(datetime.now(), ': Saved: ', medSubt_fil)

hdul[1].data = median_img
hdul[1].writeto(median_fil, overwrite=True)
print(datetime.now(), ': Saved: ', median_fil)

hdul.close()




