import matplotlib.pyplot as plt
import os
import numpy as np
from astropy.io import fits
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table, Column, join, join_skycoord
import yaml
import sys
from glob import glob
########################################################## 
def main():
    paramnames = glob('photometry_params_v*.yml')
    if len(paramnames) > 1: 
        sys.exit('error, more than one photometry_params_v*.yml')
    if len(paramnames) == 0: 
        sys.exit('no photometry_params_v*.yml file not found')

    with open(paramnames[0]) as f:
    # use safe_load instead load
        dataMap = yaml.safe_load(f)


    basedir = dataMap['basedir']
    workdir = basedir + 'photometry/sextractor/'
    convdir = basedir + 'photometry/PSF_matching/'
    stackdir  = basedir+'junk/current_best/'

    if not os.path.isdir(workdir): os.mkdir(workdir)
    os.chdir(workdir)

    qso_name = basedir.split('/')[-2]
    catname = qso_name+'_photcat_v'+str(dataMap['ver'])

    #copy config files
    os.system('cp '+basedir+'photometry/param_files/photcat/* .')

    paramfile = 'config.v4'
    det_band = dataMap['det_band']


    bands = dataMap['nrc_stacks']+dataMap['hst_stacks']
    #first make image files
    for opt in bands:
        if opt['inst'] == 'nrc': continue
        hdul = fits.open(opt['conv_name'])

        #make band images
        if opt['detection'] == True:
            masked_data = np.copy(hdul['SCI'].data)
            masked_data[(hdul['ERR'].data == 0)] = np.nan 
            fits.writeto(workdir+'detection_image_'+opt['FILTER']+'.fits',masked_data, hdul['SCI'].header, overwrite=True)
            fits.writeto(workdir+'error.fits',                           hdul['ERR'].data, hdul['SCI'].header, overwrite=True)
        else:
            if opt['inst'] == 'nrc': 
                ext = hdul['SCI']
            else:
                ext = hdul[0]
                fits.writeto(workdir+'image_'+opt['FILTER']+'_conv.fits', ext.data, ext.header, overwrite=True)


    #first make image files
    for opt in bands:
        
        #make band images
        if opt['detection'] == True:
            #run sextractor
            catname_det = workdir + 'sourcecat_det_'+opt['FILTER']+'.fits'
            if opt['inst'] == 'nrc': continue
            oargs = '-CATALOG_NAME %s' % (catname_det) +   ' -PHOT_AUTOPARAMS 2.5,3.5 -MAG_ZEROPOINT %f' % (opt['zpt'])
            cmd = 'source-extractor detection_image_'+det_band+'.fits detection_image_'+opt['FILTER']+'.fits  -c '+paramfile+' '+oargs
            #os.system(cmd)
            
            catname_phot = workdir + 'photcat_dual_'+opt['FILTER']+'.fits'
            oargs = '-CATALOG_NAME %s' % (catname_phot) +' -MAG_ZEROPOINT %f' % (opt['zpt'])
            bkgsave = ' -CHECKIMAGE_NAME  bkg_'+opt['FILTER']+'.fits '
            cmd = 'source-extractor detection_image_'+det_band+'.fits detection_image_'+opt['FILTER']+'.fits -c '+paramfile+' '+oargs+bkgsave
            #os.system(cmd)
            #print(cmd)

        else:
            catname_phot = workdir + 'photcat_dual_'+opt['FILTER']+'.fits'
            oargs = '-CATALOG_NAME %s' % (catname_phot) +' -MAG_ZEROPOINT %f' % (opt['zpt'])
            bkgsave = ' -CHECKIMAGE_NAME  bkg_'+opt['FILTER']+'.fits '
            cmd = 'source-extractor detection_image_'+det_band+'.fits image_'+opt['FILTER']+'_conv.fits -c '+paramfile+' '+oargs + bkgsave
            #os.system(cmd)
            #print(cmd)

    #load and merge cats
    source_cat_tab = Table.read(catname_det)
    for col in source_cat_tab.columns[1:]:  source_cat_tab[col].name = source_cat_tab[col].name + '_det'

    for opt in bands:
        catname_phot = workdir + 'photcat_dual_'+opt['FILTER']+'.fits'
        phot_cat = Table.read(catname_phot)
        #rename cols
        for col in phot_cat.columns[1:]:  phot_cat[col].name = phot_cat[col].name + '_' + opt['FILTER']

        #merge
        if opt['detection']:
            join_temp = join(source_cat_tab, phot_cat, keys='NUMBER', join_type='left')
        else:
            join_temp = join(join_temp,      phot_cat, keys='NUMBER', join_type='left')

    #add apcor and flux cols
    master_cat = join_temp

    master_cat.add_column(Column((master_cat['MAG_AUTO_det']-master_cat['MAG_AUTO_'+det_band])-0.090073 , name='ap_corr'))
    master_cat.add_column(Column(10.0**(-0.4*master_cat['ap_corr']) , name='ap_flux_fact'))

    #select good cols
    bad_mag = (master_cat['MAG_AUTO_det'] == 99.)
    bad_nrc = (master_cat['MAG_AUTO_det'] == 99.)
    keep_cols = ['NUMBER', 'ALPHA_J2000_det', 'DELTA_J2000_det', 'X_IMAGE_det', 'Y_IMAGE_det']
    for opt in bands:
        #add mag apcor
        obsflux_to_njy = 10.**(9.0+0.4*8.9-0.4*opt['zpt']) 
        master_cat.add_column(Column((master_cat['MAG_AUTO_'+opt['FILTER']]+master_cat['ap_corr']), name='MAG_AUTO_'+opt['FILTER']+'_apcor'))  
        #add flux for EAZY
        master_cat.add_column(Column(master_cat['FLUX_AUTO_'+opt['FILTER']]*master_cat['ap_flux_fact']*obsflux_to_njy , name='fnu_'+opt['FILTER']+'_AUTO_apcor')) 
        #add ers
        master_cat.add_column(Column(master_cat['FLUXERR_AUTO_'+opt['FILTER']]*master_cat['ap_flux_fact']*obsflux_to_njy , name='enu_'+opt['FILTER']+'_AUTO_apcor'))  #add to AUTO mags
        
        #flag for good mags
        bad_mag = (bad_mag)|(master_cat['MAG_AUTO_'+opt['FILTER']] == 99.)
        if opt['inst'] == 'nrc':
            bad_nrc = (bad_nrc)|(master_cat['MAG_AUTO_'+opt['FILTER']] == 99.)
        keep_cols = keep_cols + ['MAG_AUTO_'+opt['FILTER']+'_apcor', 'MAGERR_AUTO_'+opt['FILTER'],'fnu_'+opt['FILTER']+'_AUTO_apcor', 'enu_'+opt['FILTER']+'_AUTO_apcor']


    #add bad obj col
    master_cat.add_column(Column(~bad_mag ,  name='good_mag')) 
    master_cat.add_column(Column(~bad_nrc ,  name='good_nrc')) 
    keep_cols = keep_cols + ['good_mag', 'good_nrc']
    #write cat
    master_cat.write(workdir + catname + '_full.fits', overwrite=True)
    ##trim down table and add flux cols
    #master_cat.keep_columns(keep_cols)
    #master_cat.write(workdir + catname + '_short_EAZY.fits', overwrite=True)


##########################################################  
if __name__ == "__main__":
    main()

