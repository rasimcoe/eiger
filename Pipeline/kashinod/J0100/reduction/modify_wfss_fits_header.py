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
pointings=[1,2,3,4,5,6,7,8,9,10,11,12]

inp_fils = []
for pt in pointings:
    #tmp = sorted(glob.glob(inp_dir+'jw01243001001_0[2,4]101_'+str(pt).zfill(5)+'_nrc*long_rate.fits'))
    #tmp = sorted(glob.glob(inp_dir+'jw01243001002_0[2,4]101_'+str(pt).zfill(5)+'_nrc*long_rate.fits'))
    tmp = sorted(glob.glob(inp_dir+'jw01243001003_0[2,4]101_'+str(pt).zfill(5)+'_nrc*long_rate.fits'))
    #tmp = sorted(glob.glob(inp_dir+'jw01243001004_0[2,4]101_'+str(pt).zfill(5)+'_nrc*long_rate.fits'))
    

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
