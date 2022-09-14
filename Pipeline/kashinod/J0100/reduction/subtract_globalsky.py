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


# This script is for processing Image2 cal.fits files
# Usage: python subtract_globalsky inp_fil out_dir globalsky_fil

# Input fits file
inp_fil = sys.argv[1]
fits_name = inp_fil.split("/")[-1]

# Output directory
out_dir = sys.argv[2]
if os.path.isdir(out_dir) == False:
    print('Make directory: ', out_dir, flush=True)
    os.mkdir(out_dir)
else:
    print('Output directory already exits: ', out_dir, flush=True)

# Global sky image
global_sky_fil = sys.argv[3]
if os.path.isfile(global_sky_fil) == False:
    raise ValueError('No global sky file exists: '+ global_sky_fil)


# Output file name
# The median-subtracted file has the same name as the input 
# for subsequence reduction processes
gskySubt_fil = os.path.join(out_dir, fits_name)
plt_fil = gskySubt_fil.replace(".fits",".pdf")


print('Input fits file : ', inp_fil, flush=True)
print('Global sky fits file : ', global_sky_fil, flush=True)
print('Output globalsky-subtracted fits file: ', gskySubt_fil, flush=True)

if inp_fil == gskySubt_fil:
    raise ValueError('Output file name must be different from the input.  Change the different output directory (', fits_name, ')', flush=True)


# Open global sky fits file
global_sky = fits.getdata(global_sky_fil)

# Open fits file
# EXT 1 SCI
# EXT 2 ERR
# EXT 3 DQ
# EXT 4 VAR_POISSON
# EXT 5 VAR_RNOISE
hdul = fits.open(inp_fil)
sci_img = np.copy(hdul[1].data)

hdul[1].data = sci_img - global_sky
hdul.writeto(gskySubt_fil, overwrite=True)
print(datetime.now(), ': Saved: ', gskySubt_fil)
hdul.close()


print(plt_fil, flush=True)


# Plot
sci_ptiles = np.nanpercentile(sci_img[:,:], [0.15, 1., 2.5, 16., 50., 84., 97.5, 99., 99.85])
print('SCI Ptiles: ', sci_ptiles)
vmin=sci_ptiles[2]
vmax=sci_ptiles[-3]

fig = plt.figure(figsize=(20,20))
ax = fig.subplots(2)
ax[0].set_title('before subtraction')
im1=ax[0].imshow(sci_img, vmin=vmin, vmax=vmax, origin="lower")

ax[1].set_title('after subtraction')
im2=ax[1].imshow(sci_img - global_sky, 
                 vmin=vmin-np.nanmedian(global_sky), 
                 vmax=vmax-np.nanmedian(global_sky), 
                 origin="lower")

fig.colorbar(im1, cax=mpl_toolkits.axes_grid1.make_axes_locatable(ax[0]).append_axes('right','5%',pad='3%'))
fig.colorbar(im2, cax=mpl_toolkits.axes_grid1.make_axes_locatable(ax[1]).append_axes('right','5%',pad='3%'))
fig.savefig(plt_fil)


