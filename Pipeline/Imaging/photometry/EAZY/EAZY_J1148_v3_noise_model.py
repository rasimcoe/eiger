import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import fits
import eazy 
import warnings
from astropy.utils.exceptions import AstropyWarning
from astropy.table import Table      
import eazy.hdf5                                                                                        
# Symlink templates & filters from the eazy-code repository
np.seterr(all='ignore')
warnings.simplefilter('ignore', category=AstropyWarning)

basedir = '/scratch/mruari/EIGER/imaging/J1148+5251/'                                       #change this
workdir = basedir + 'photometry/EAZY/'
os.chdir(workdir)
if not os.path.exists('templates'):
    eazy.symlink_eazy_inputs() 

basename = workdir+'J1148_photzcat_v2_EAZY_'                                                #change this

cat = Table.read(basedir+"photometry/catalogs/J1148+5251_photcat_v2_noisemodel_short.fits") #change this
#cat = cat[(cat['good_mag'])] 

params = {}
# Galactic extinction
#params['MW_EBV'] = 0.00
#params['CAT_HAS_EXTCORR'] = True
params['Z_STEP'] = 0.01
params['Z_MIN'] = 0.05
params['Z_MAX'] = 10.
#params['SYS_ERR'] = 0.03

params['APPLY_PRIOR'] = False                   #does nothing
params['PRIOR_ABZP'] = 31.4         #nJy         #microJY 23.9
params['PRIOR_FILTER'] = 366 # K
params['PRIOR_FILE'] = 'templates/prior_K_TAO.dat'

params['TEMPLATES_FILE'] = 'templates/fsps_full/tweak_fsps_QSF_12_v3.param'
params['TEMPLATE_COMBOS'] = 1
params['FIX_ZSPEC'] = False
params['IGM_SCALE_TAU'] = 1.0
params['N_MIN_COLORS'] = 3
params['NOT_OBS_THRESHOLD'] = -90
params['CAT_HAS_EXTCORR'] = False
params['MW_EBV'] = 0.0273                                                           #from Planck Change dust extinction for different field

#inputs
cat_name = basename + '_input.csv'
params['CATALOG_FILE'] = cat_name
params['CATALOG_FORMAT'] = 'csv'

#output
params['MAIN_OUTPUT_FILE'] = basename + '_output'

# Get ttable
#make cat with new needed cols
cat_fmt = Table([cat['NUMBER'], 
        cat['F606W_AUTO_fnu'].data, cat['F606W_AUTO_enu'].data,                     #Change these lines for different HST bands
        cat['F775W_AUTO_fnu'].data, cat['F775W_AUTO_enu'].data,                     #Change these lines for different HST bands
        cat['F850LP_AUTO_fnu'].data,cat['F850LP_AUTO_enu'].data,                    #Change these lines for different HST bands
        cat['F115W_AUTO_fnu'].data, cat['F115W_AUTO_enu'].data,
        cat['F200W_AUTO_fnu'].data, cat['F200W_AUTO_enu'].data,
        cat['F356W_AUTO_fnu'].data, cat['F356W_AUTO_enu'].data,
        np.zeros(len(cat), dtype=float)], 
        names=('id',
             'HST_F606W_fnu','HST_F606W_enu',                                       #Change these lines for different HST bands
             'HST_F775W_fnu','HST_F775W_enu',                                       #Change these lines for different HST bands
             'HST_F850LP_fnu','HST_F850LP_enu',                                     #Change these lines for different HST bands
             'NRC_F115W_fnu','NRC_F115W_enu',
             'NRC_F200W_fnu','NRC_F200W_enu',
             'NRC_F356W_fnu','NRC_F356W_enu',
             'z_spec'))

cat_fmt.write(cat_name, format='csv', overwrite=True)


# translate file
curnames = [ 'HST_F606W_fnu','HST_F606W_enu',                                       #Change these lines for different HST bands
             'HST_F775W_fnu','HST_F775W_enu',                                       #Change these lines for different HST bands
             'HST_F850LP_fnu','HST_F850LP_enu',                                     #Change these lines for different HST bands
             'NRC_F115W_fnu','NRC_F115W_enu',
             'NRC_F200W_fnu','NRC_F200W_enu',
             'NRC_F356W_fnu','NRC_F356W_enu']

transnames =['F236','E236',                                       #Change these lines for different HST bands       #236    173 hst/ACS_update_sep07/wfc_f606w_t81.dat obs_AFTER_7-4-06+rebin-5A lambda_c= 5.9211e+03 AB-Vega= 0.083 w95=2224.8
             'F238','E238',                                       #Change these lines for different HST bands       #238     86 hst/ACS_update_sep07/wfc_f775w_t81.dat obs_AFTER_7-4-06+rebin-5A lambda_c= 7.6924e+03 AB-Vega= 0.384 w95=1490.9
             'F240','E240',                                       #Change these lines for different HST bands       #240    102 hst/ACS_update_sep07/wfc_f850lp_t81.dat obs_AFTER_7-4-06+rebin-5A lambda_c= 9.0331e+03 AB-Vega= 0.517 w95=2091.9 ,
             'F364','E364',
             'F366','E366',
             'F376','E376']
        

temp_cat = Table([curnames, transnames],names=('column','trans'))
trans_name = basename + '_translate.csv'
temp_cat.write(trans_name, format='csv', overwrite=True)


self = eazy.photoz.PhotoZ(param_file=None, translate_file=trans_name, zeropoint_file=None, 
                          params=params, load_prior=True, load_products=False)

#apply SFH constraint
self.tempfilt.apply_SFH_constraint()

# #Iterative zeropoint corrections
# NITER = 3
# NBIN = np.minimum(self.NOBJ//100, 180)
# self.param.params['VERBOSITY'] = 1.
# for iter in range(NITER):
#     print('Iteration: ', iter)
    
#     sn = self.fnu/self.efnu
#     clip = (sn > 1).sum(axis=1) > 2 # Generally make this higher to ensure reasonable fits
#     self.iterate_zp_templates(idx=self.idx[clip], update_templates=False, 
#                               update_zeropoints=True, iter=iter, n_proc=8, 
#                               save_templates=False, error_residuals=False, 
#                               NBIN=NBIN, get_spatial_offset=False)

# # Turn off error corrections derived above
# self.set_sys_err(positive=True)

# Full catalog
sample = np.isfinite(self.ZSPEC)
#open arrays for loop
ids = self.idx
z_spec = self.ZSPEC

# fit_parallel renamed to fit_catalog 14 May 2021
self.fit_catalog(self.idx[sample], n_proc=8, prior=False, beta_prior=False)                #, prior=False, beta_prior=False

zout, hdu = self.standard_output(simple=False, 
                                 rf_pad_width=0.5, rf_max_err=2, 
                                 prior=False, beta_prior=False, 
                                 absmag_filters=[], 
                                 extra_rf_filters=[])

eazy.hdf5.write_hdf5(self, h5file=self.param['MAIN_OUTPUT_FILE'] + '.h5')

self.show_fit(9699-1, xlim=[0.2, 10], show_components=True, show_prior=True, logpz=True, zr=[0,10])