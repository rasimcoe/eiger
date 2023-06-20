import matplotlib.pyplot as plt
import os
import sys
import numpy as np
from astropy.io import fits
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table, Column, join, join_skycoord
from photutils import CircularAperture, aperture_photometry, CircularAnnulus, EllipticalAperture, ApertureStats, EllipticalAnnulus
from glob import glob
import yaml
###############################################################
def measure_mean_WHT(WHT, catalog, FILTER):
    #good sources
    obj_wht = np.zeros(len(catalog), dtype=float)
    for io, obj in enumerate(catalog):
        posn = np.array([obj['X_IMAGE_det']-1, obj['Y_IMAGE_det']-1])
        if (obj['A_IMAGE_det'] > 0.0)&(obj['B_IMAGE_det'] > 0.0)&(obj['KRON_RADIUS_det'] > 0.0):
            r_kron = obj['KRON_RADIUS_det'] * 1.2/2.5
            a_image = (obj['A_IMAGE_det']*r_kron)
            b_image = (obj['B_IMAGE_det']*r_kron)
            theta_image = obj['THETA_IMAGE_det'] * u.deg
            apers = EllipticalAperture(posn.T, a_image, b_image, theta=theta_image)
            
            phot_table = aperture_photometry(WHT, apers)
            obj_wht[io] = phot_table['aperture_sum'].data / apers.area
        else:
            #just closest pixel
            obj_wht[io] = WHT[np.round(posn).astype(int)[1],np.round(posn).astype(int)[0]]
    scaled_wht = obj_wht/np.median(WHT[(WHT != 0)])
    return scaled_wht
###############################################################
def main():
    paramnames = glob('photometry_params_v*.yml')
    if len(paramnames) > 1: 
        sys.exit('error, more than one photometry_params_v*.yml')
    if len(paramnames) == 0: 
        sys.exit('no photometry_params_v*.yml file not found')

    with open(paramnames[0]) as f:
    # use safe_load instead load
        dataMap = yaml.safe_load(f)


    nrc_stacks = dataMap['nrc_stacks']
    hst_stacks = dataMap['hst_stacks']

    nrc_stacks = [nrc_stacks[0],nrc_stacks[2],nrc_stacks[1]]
    hst_stacks = [hst_stacks[0]]

    basedir = dataMap['basedir']
    workdir = basedir + 'photometry/sextractor/'
    convdir = basedir + 'photometry/PSF_matching/'
    noisedir = basedir + 'photometry/noise_model/'
    stackdir  = basedir+'junk/current_best/'
    catdir = basedir + 'photometry/catalogs/'
    if not os.path.isdir(catdir): os.mkdir(catdir)
    os.chdir(workdir)



    qso_name = basedir.split('/')[-2]
    catname = qso_name+'_photcat_v'+str(dataMap['ver'])+'_full.fits'
    catalog = Table.read(catname)

    save_name = catdir + qso_name+'_photcat_v'+str(dataMap['ver'])+'_noisemodel_full.fits'
    short_name = catdir + qso_name+'_photcat_v'+str(dataMap['ver'])+'_noisemodel_short.fits'

    noise_model = np.load(noisedir + 'noise_model_params_v%s.npz'%(dataMap['ver']))

    ap_diam = ([5,7,11,17])

    det_cols = [
                 'NUMBER',
                 'ID_PARENT_det',
                 'ALPHA_J2000_det',
                 'DELTA_J2000_det',
                 'X_IMAGE_det',
                 'Y_IMAGE_det',
                 'A_IMAGE_det',
                 'B_IMAGE_det',
                 'THETA_IMAGE_det',
                 'ELONGATION_det',
                 'ELLIPTICITY_det',
                 'FWHM_WORLD_det',
                 'CLASS_STAR_det',
                 'BACKGROUND_det',
                 'ap_corr',
                 'good_mag',
                 'good_nrc']

    mauto_cols = []
    fauto_cols = []
    map_cols   = []
    fap_cols   = []


    for istack, stack in enumerate(nrc_stacks):
        FILTER = stack['FILTER']
        #add column
        catalog.add_column(np.zeros(len(catalog), dtype=float), name='enu_'+FILTER+'_aper_model')

        hdul = fits.open(stack['conv_name'])

        #calc obj wht radios
        cat_scale_wht = measure_mean_WHT(hdul['WHT'].data, catalog, FILTER)

        #calc area
        r_kron = catalog['KRON_RADIUS_det'] * 1.2/2.5
        auto_area = (np.pi *
        (catalog['A_IMAGE_det']*r_kron).value *
        (catalog['B_IMAGE_det']*r_kron).value)

        #calc noise
        isave = np.where(noise_model['filters'] == FILTER)[0][0]
        noise_uni = noise_model['popt'][isave,0] * auto_area**(noise_model['popt'][isave,1])
        #scale to wht
        noise_scale = noise_uni / np.sqrt(cat_scale_wht)

        #add to catalog
        catalog[FILTER+'_AUTO_enu'] = noise_scale * (catalog['fnu_'+FILTER+'_AUTO_apcor']/catalog['FLUX_AUTO_'+FILTER])
        catalog[FILTER+'_AUTO_fnu'] = catalog['fnu_'+FILTER+'_AUTO_apcor']
        catalog[FILTER+'_AUTO_mag'] = catalog['MAG_AUTO_'+FILTER+'_apcor']

        mauto_cols = mauto_cols + [FILTER+'_AUTO_mag']
        fauto_cols = fauto_cols + [FILTER+'_AUTO_fnu', FILTER+'_AUTO_enu']


        obsflux_to_njy = 10.**(9.0+0.4*8.9-0.4*stack['zpt']) 
        for iap in range(len(ap_diam)):
            ap_area = np.pi * ap_diam[iap]**2
            noise_uni = noise_model['popt'][isave,0] * ap_area**(noise_model['popt'][isave,1])
            noise_scale = noise_uni / np.sqrt(cat_scale_wht)
            catalog[FILTER+'_AP%i_enu'%(iap+1)] = noise_scale * obsflux_to_njy
            catalog[FILTER+'_AP%i_fnu'%(iap+1)] = catalog['FLUX_APER_'+FILTER][:,iap] * obsflux_to_njy
            catalog[FILTER+'_AP%i_mag'%(iap+1)] = catalog['MAG_APER_'+FILTER][:,iap]
            map_cols = map_cols + [FILTER+'_AP%i_mag'%(iap+1)]
            fap_cols = fap_cols + [FILTER+'_AP%i_fnu'%(iap+1), FILTER+'_AP%i_enu'%(iap+1)]



    for istack, stack in enumerate(hst_stacks):
        FILTER = stack['FILTER']

        #add column
        catalog.add_column(np.zeros(len(catalog), dtype=float), name='enu_'+FILTER+'_aper_model')

        hdul = fits.open(stack['conv_name'])

        #calc area
        r_kron = catalog['KRON_RADIUS_det'] * 1.2/2.5
        auto_area = (np.pi *
        (catalog['A_IMAGE_det']*r_kron).value *
        (catalog['B_IMAGE_det']*r_kron).value)

        #calc noise
        isave = np.where(noise_model['filters'] == FILTER)[0][0]
        noise_uni = noise_model['popt'][isave,0] * auto_area**(noise_model['popt'][isave,1])

        #add to catalog, scale to nJy
        #add to catalog
        catalog[FILTER+'_AUTO_enu'] = noise_uni * (catalog['fnu_'+FILTER+'_AUTO_apcor']/catalog['FLUX_AUTO_'+FILTER])
        catalog[FILTER+'_AUTO_fnu'] = catalog['fnu_'+FILTER+'_AUTO_apcor']
        catalog[FILTER+'_AUTO_mag'] = catalog['MAG_AUTO_'+FILTER+'_apcor']

        mauto_cols = mauto_cols + [FILTER+'_AUTO_mag']
        fauto_cols = fauto_cols + [FILTER+'_AUTO_fnu', FILTER+'_AUTO_enu']

        obsflux_to_njy = 10.**(9.0+0.4*8.9-0.4*stack['zpt']) 
        for iap in range(len(ap_diam)):
            ap_area = np.pi * ap_diam[iap]**2
            noise_uni = noise_model['popt'][isave,0] * ap_area**(noise_model['popt'][isave,1])
            catalog[FILTER+'_AP%i_enu'%(iap+1)] = noise_uni * obsflux_to_njy
            catalog[FILTER+'_AP%i_fnu'%(iap+1)] = catalog['FLUX_APER_'+FILTER][:,iap] * obsflux_to_njy
            catalog[FILTER+'_AP%i_mag'%(iap+1)] = catalog['MAG_APER_'+FILTER][:,iap]
            map_cols = map_cols + [FILTER+'_AP%i_mag'%(iap+1)]
            fap_cols = fap_cols + [FILTER+'_AP%i_fnu'%(iap+1), FILTER+'_AP%i_enu'%(iap+1)]

    keep_cols = det_cols + mauto_cols + fauto_cols + map_cols + fap_cols

    #flag bad fluxes with large negative values
    filter_list = [test['FILTER'] for  test in nrc_stacks+hst_stacks]
    for FILTER in filter_list:
        bad_flux = (catalog['fnu_'+FILTER+'_AUTO_apcor'] == 0.)
        catalog[FILTER+'_AUTO_enu'][(bad_flux)]   = -99
        catalog[FILTER+'_AUTO_fnu'][(bad_flux)]   = -99
        for iap in range(len(ap_diam)):
            catalog[FILTER+'_AP%i_enu'%(iap+1)][(bad_flux)]   = -99
            catalog[FILTER+'_AP%i_fnu'%(iap+1)][(bad_flux)]    = -99

    catalog.write(save_name, overwrite=True)
    catalog = catalog[keep_cols]
    catalog.write(short_name, overwrite=True)

##########################################################  
if __name__ == "__main__":
    main()
