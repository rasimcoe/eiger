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

import pyregion
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS, utils


# Usage: python mask_objects_img2cal.py inp_fil out_dir output_mask_fil region_fil radius_factor
print(datetime.now(), '- Start.')
# Input fits file
inp_fil = sys.argv[1]
fits_name = inp_fil.split("/")[-1]
print('Input fits file: ', fits_name, flush=True)

# Output directory
out_dir = sys.argv[2]
if os.path.isdir(out_dir) == False:
    print('Make directory: ', out_dir, flush=True)
    os.mkdir(out_dir)
else:
    print('Output directory already exits: ', out_dir, flush=True)


# Output mask file name
mask_fil = fits_name.replace("_cal.fits", "_mask.fits")
mask_fil = os.path.join(out_dir, mask_fil)
print('Output mask fits file: ', mask_fil, flush=True)

# Region file
reg_fil = sys.argv[3]
if os.path.isfile(reg_fil) == False:
    raise ValueError('No region file exists: '+ reg_fil)

print('Read region file: ', reg_fil, flush=True)
print(datetime.now(), flush=True)
print('This takes a couple of minutes. Wait.', flush=True)
region_list = pyregion.open(reg_fil)
print(datetime.now(), flush=True)
print('Number of regions: ', len(region_list), flush=True)

# Increase the radius by a factor of RAD_FACT
if len(sys.argv) == 5:
    rad_fact = float(sys.argv[4])
    print('Increase the radii of elliptical regions by ', rad_fact, flush=True)
    region_list_orig = region_list.copy()
    for i in range(len(region_list)):
        if region_list[i].name == 'ellipse':
            #print('--------', flush=True)
            #print(region_list[i].coord_list, flush=True)
            region_list[i].coord_list[2] = region_list_orig[i].coord_list[2] * rad_fact
            region_list[i].coord_list[3] = region_list_orig[i].coord_list[3] * rad_fact
            #print(region_list[i].coord_list, flush=True)

# Open fits file
# EXT 1 SCI
# EXT 2 ERR
# EXT 3 DQ
# EXT 4 VAR_POISSON
# EXT 5 VAR_RNOISE

print('Read '+inp_fil, flush=True)
hdul = fits.open(inp_fil)
header = hdul[1].header

print(datetime.now(), ' -- Get mask.  This takes some time. Wait.', flush=True)
mask = region_list.get_mask(hdu=hdul[1])
print(datetime.now(), ' -- Get_mask done.', flush=True)

hdul.close()
mask_ = np.zeros(mask.shape, dtype="int32")
mask_[mask]=1

hdu = fits.PrimaryHDU(mask_, header)
hdul = fits.HDUList([hdu])
hdul.writeto(mask_fil, overwrite=True)
print('Saved: '+mask_fil, flush=True)

hdul.close()
print(datetime.now(), '- Finished.')

