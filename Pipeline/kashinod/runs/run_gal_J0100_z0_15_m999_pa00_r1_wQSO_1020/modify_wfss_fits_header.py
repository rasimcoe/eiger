import glob
import os
import numpy as np
from astropy.io import fits

inp_dir='./calibrated_det1/'
out_dir='./calibrated_det1_wfss_hdr_corr/'

if not os.path.isdir(out_dir):
    print('make directory: ', out_dir)
    os.mkdir(out_dir)

# WFSS pointings
pointings = [1,2,3,4,5,6,7,8,9,10,11,12,
             14,15,16,17,18,19,20,21,22,23,24,25,
             29,30,31,32,33,34,35,36,37,38,39,40,
             42,43,44,45,46,47,48,49,50,51,52,53,
             57,58,59,60,61,62,63,64,65,66,67,68,
             70,71,72,73,74,75,76,77,78,79,80,81,
             85,86,87,88,89,90,91,92,93,94,95,96,
             98,99,100,101,102,103,104,105,106,107,108,109]

inp_fils = []
for pt in pointings:
    tmp = sorted(glob.glob(inp_dir+'jw0124300100?_01101_'+str(pt).zfill(5)+'_nrc*5_rate.fits'))
    inp_fils.extend(tmp)

for inp_fil in inp_fils:

    print('-------------------------------')
    print('Input file: ', inp_fil)

#inp_fil = os.path.join(inp_dir, fits_name)
    fits_name=inp_fil.split("/")[-1]
    out_fil = os.path.join(out_dir, fits_name)
    print('Output file: ', out_fil)
    if inp_fil==out_fil:
        raise ValueError
    
    hdul = fits.open(inp_fil)
    print('PUPIL   : ', hdul[0].header['pupil'])
    print('EXP_TYPE: ', hdul[0].header['exp_type'])

    print('Modigy header...')
    hdul[0].header['pupil'] = 'CLEAR'
    hdul[0].header['exp_type'] = 'NRC_IMAGE'
    print('PUPIL   : ', hdul[0].header['pupil'])
    print('EXP_TYPE: ', hdul[0].header['exp_type'])

    hdul.writeto(out_fil, overwrite=True)
    hdul.close()
