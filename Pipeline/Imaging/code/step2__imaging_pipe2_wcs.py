import shutil
from glob import glob
import os
os.environ["CRDS_DATA"] = "/scratch/mruari/EIGER/cache/crds_cache"
os.environ["CRDS_PATH"] = "/scratch/mruari/EIGER/cache/crds_cache"
os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"
import sys
import logging
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.nddata import NDData
from astroquery.mast import Observations
from photutils import detect_threshold, DAOStarFinder, IRAFStarFinder
from stwcs.wcsutil import HSTWCS
#from drizzlepac import updatehdr
from astropy.table import Table, Column
from tweakwcs import fit_wcs, align_wcs, FITSWCS, TPMatch, WCSImageCatalog
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
from jwst.resample.resample_utils import build_mask
from jwst import datamodels
from jwst.tweakreg.tweakreg_step import make_tweakreg_catalog, _common_name
from jwst.datamodels import ImageModel
from tweakwcs import JWSTgWCS
from jwst.pipeline import Image3Pipeline
from jwst.tweakreg import astrometric_utils as amutils
from photutils.psf import extract_stars
from multiprocessing import Pool
import yaml
################################################
def make_star_cat(ref_img, FILTER):
    
    mask = (ref_img['WHT'].data == 0)|(ref_img['SCI'].data == 0)

    mean, median, std = sigma_clipped_stats(ref_img['SCI'].data, sigma=3.0, mask=mask)  

    daofind = DAOStarFinder(fwhm=1.8, threshold=200.*std,
                            sharplo=0.5, sharphi=1.0, roundlo=-1.0,
                            roundhi=1.0, brightest=None,
                            peakmax=100)

    sources = daofind(ref_img['SCI'].data - median, mask=mask)
    sources.rename_column('xcentroid','x')
    sources.rename_column('ycentroid','y')

    #convert to ra dec
    wcs = WCS(ref_img['SCI'].header)
    ra, dec = wcs.all_pix2world(sources['x'], sources['y'],0)
    sources.add_column(Column(ra, name='RA'))
    sources.add_column(Column(dec, name='DEC'))

    nddata = NDData(ref_img['SCI'].data)  
    all_stars = extract_stars(nddata, sources, size=31) 
    saturated = np.any(np.array(all_stars.data) == 0, axis=(1,2))
    clean_sources = sources[(~saturated)]


    #remove close dups
    #keep_sources = isolated(sources, 2.0)
    columns = ['RA', 'DEC', 'id', 'flux']
    catalog = clean_sources[columns]
    clean_sources.write('star_cat_'+FILTER+'.fits', overwrite=True)
    #sources.write('ref_cat_F356W_v4.txt', format='ascii.commented_header', overwrite=True)
    return clean_sources
################################################
def isolated(cat, minsep):
    coord_sex = SkyCoord(ra=cat['RA'], dec=cat['DEC'], unit='deg')
    #idx, d2d, d3d = coord_sex.match_to_catalog_sky(coord_sex)
    #max_sep = minsep * u.arcsec
    #sep_constraint = d2d < max_sep
    #sex_matches = sex_cat[idx[(sep_constraint)]]
    idx1, idx2, d2d, d3d = coord_sex.search_around_sky(coord_sex, minsep * u.arcsec)
    n_matches, bins = np.histogram(idx1, bins=np.arange(len(cat)+1))
    return cat[(n_matches == 1)]
################################################
def make_refcat(ref_img, FILTER):
    
    mask = (ref_img['WHT'].data == 0)

    mean, median, std = sigma_clipped_stats(ref_img['SCI'].data, sigma=3.0, mask=mask)  

    daofind = DAOStarFinder(fwhm=1.8, threshold=20.*std,
                            sharplo=0.5, sharphi=1.0, roundlo=-1.0,
                            roundhi=1.0, brightest=None,
                            peakmax=None)

    # daofind = IRAFStarFinder(fwhm=3.8, threshold=20.*std,
    #                         sharplo=0.5, sharphi=1.0, roundlo=-1.0,
    #                         roundhi=1.0, brightest=None,
    #                         peakmax=None, minsep_fwhm=10)

    # Mask the non-imaging area (e.g. MIRI)


    sources = daofind(ref_img['SCI'].data - median, mask=mask)

    #convert to ra dec
    wcs = WCS(ref_img['SCI'].header)
    ra, dec = wcs.all_pix2world(sources['xcentroid'], sources['ycentroid'],0)
    sources.add_column(Column(ra, name='RA'))
    sources.add_column(Column(dec, name='DEC'))

    sources.rename_column('xcentroid','x')
    sources.rename_column('ycentroid','y')

    #remove close dups
    keep_sources = isolated(sources, 2.0)
    columns = ['RA', 'DEC', 'id', 'flux']
    catalog = keep_sources[columns]
    keep_sources.write('ref_cat_'+FILTER+'.fits', overwrite=True)
    #sources.write('ref_cat_F356W_v4.txt', format='ascii.commented_header', overwrite=True)
    return keep_sources
################################################
def update_fits_hdr(dm):
    crpix = [dm.meta.wcsinfo.crpix1, dm.meta.wcsinfo.crpix2]
    crpix = None if None in crpix else crpix
    fit_sip_hdr = dm.meta.wcs.to_fits_sip(
    max_pix_error=1e-4, # adjust accuracy and degree to your needs
    degree=range(1, 4),
    max_inv_pix_error=1e-4,
    inv_degree=range(1, 4),
    npoints=400,
    crpix=crpix
    )

    # update meta.wcs_info with fit keywords except for naxis*
    del fit_sip_hdr['naxis*']

    # maintain convention of lowercase keys
    fit_sip_hdr = {k.lower(): v for k, v in fit_sip_hdr.items()}

    # delete naxis, cdelt, pc from wcsinfo
    rm_keys = ['naxis', 'cdelt1', 'cdelt2', 'pc1_1', 'pc1_2', 'pc2_1', 'pc2_2']
    for key in rm_keys:
        if key in dm.meta.wcsinfo.instance:
            del dm.meta.wcsinfo.instance[key]

    # update meta.wcs_info with fit keywords
    dm.meta.wcsinfo.instance.update(fit_sip_hdr) 
    return dm
###############################################################################
def apply_transform(fix_cal, temp_dir, aligned_imwcs, ref_wcs=None):
    for filename in fix_cal:
        dm = ImageModel(filename)
        # create JWST corrector object using GWCS from the data model:
        tpwcs_corrector = JWSTgWCS(dm.meta.wcs, dm.meta.wcsinfo.instance)
        # apply some known correction (typically use fit_wcs or align_wcs instead
        # of directly setting corrections the way it is done below): NO TRANSPOSE on matrix
        tpwcs_corrector.set_correction(matrix=aligned_imwcs.meta['matrix'], shift=-aligned_imwcs.meta['shift'], ref_tpwcs=ref_wcs)
        # assign corrected WCS back to the data model:
        dm.meta.wcs = tpwcs_corrector.wcs
        #update fits hdr
        dm = update_fits_hdr(dm)
        # save image model to a different file but one can also overwrite existing model:
        dm.write(temp_dir + filename.split('/')[-1]) 
###############################################################################
def gaia_align(basedir, FILTER, restack=False):

    cal_dir = basedir+FILTER+'/pipe1_basic/'
    temp_dir = basedir+FILTER+'/pipe2_mywcs/'
    if not os.path.isdir(temp_dir): os.mkdir(temp_dir)
    os.chdir(temp_dir)

    #load ref image

    img_hudl = fits.open(basedir + 'F356W/pipe1_basic/stack_F356W_pipe1.fits')
    #load wcs
    img_w = HSTWCS(img_hudl, ('SCI', 1))
    img_wcs = FITSWCS(img_w)
    #make ref cat
    img_cat = make_star_cat(img_hudl, FILTER)

    #find center ra dec
    cen_coord = img_wcs.wcs.pixel_to_world(img_hudl['SCI'].data.shape[1]/2., img_hudl['SCI'].data.shape[0]/2.)

    gaia_cat = amutils.get_catalog(cen_coord.ra.value, cen_coord.dec.value, sr=0.1, catalog='GAIADR2')
    colnames = ('ra', 'dec', 'mag', 'objID')
    ref_table = gaia_cat[colnames]
    ref_table.rename_column('ra','RA')
    ref_table.rename_column('dec','DEC')
    ref_table.write('Gaia_ref_table.fits')


    match = TPMatch(searchrad=10, separation=0.5, tolerance=5, use2dhist=True)
    ridx, iidx = match(ref_table, img_cat, img_wcs)

    # Align image WCS:
    aligned_imwcs = fit_wcs(ref_table[ridx], img_cat[iidx], img_wcs, fitgeom='rshift')

    #find cal files
    fix_cal = glob(cal_dir+'jw*_crf.fits')
    #apply transform
    apply_transform(fix_cal, temp_dir, aligned_imwcs, ref_wcs=img_wcs)

    # #restack
    # if restack == True:
    #     tweak_files = glob(temp_dir+'jw*_crf.fits')
    #     stackname = 'stack_'+FILTER+'_pipe2_mywcs_rotfix.fits'
    #     restack(tweak_files, basedir, temp_dir, FILTER, stackname)
###############################################################################
def internal_align(basedir, FILTER, restack=False):
    #basedir, basedir = '/scratch/mruari/EIGER/imaging/J0100+2802/'
    cal_dir = basedir+FILTER+'/pipe1_basic/'
    temp_dir = basedir+FILTER+'/pipe2_mywcs/'
    if not os.path.isdir(temp_dir): os.mkdir(temp_dir)
    os.chdir(temp_dir)

    #load filt image
    ref_img = fits.open(basedir + 'F356W/pipe2_mywcs/stack_F356W_pipe2_mywcs.fits')
    # load wcs of img to be aligne
    ref_w = HSTWCS(ref_img, ('SCI', 1))
    ref_wcs = FITSWCS(ref_w)
    #make align cat
    ref_cat = make_refcat(ref_img, 'F356W')
    
    #load filt image
    align_img = fits.open(cal_dir + 'stack_'+FILTER+'_pipe1.fits')
    # load wcs of img to be aligne
    align_w = HSTWCS(align_img, ('SCI', 1))
    align_wcs = FITSWCS(align_w)
    #make align cat
    align_cat = make_refcat(align_img, FILTER)

    # Match sources in the catalogs:
    match = TPMatch(searchrad=10, separation=0.5, tolerance=5, use2dhist=True)
    ridx, iidx = match(ref_cat, align_cat, align_wcs)
    #solve wcs
    aligned_imwcs = fit_wcs(ref_cat[ridx], align_cat[iidx], align_wcs, fitgeom='rshift')

    #find cal files
    fix_cal = glob(cal_dir+'jw*_crf.fits')
    #apply transform
    apply_transform(fix_cal, temp_dir, aligned_imwcs, ref_wcs)

    # #restack
    # if restack == True:
    #     tweak_files = glob(temp_dir+'jw*_crf.fits')
    #     stackname = 'stack_'+FILTER+'_pipe2_mywcs_rotfix.fits' 
    #     restack(tweak_files, basedir, temp_dir, FILTER, stackname)
#################################################################################
def restack(params):
    #load wcs 
    FILTER = params['FILTER']
    basedir = params['basedir']
    save_dir = basedir+FILTER+'/pipe2_mywcs/'

    #load wcs 
    wcs_pickle = basedir + 'output_grid_wcs.npz'
    save_wcs = np.load(wcs_pickle)

    #Step 3:
    fits_files_corr = glob(save_dir + 'jw*crf.fits')

    pipe = Image3Pipeline()
    pipe.output_dir=save_dir
    pipe.tweakreg.skip=True
    pipe.outlier_detection.skip=True
    pipe.save_results=True
    pipe.skymatch.skip=True
    #fixed grid for resample
    pipe.resample.pixel_scale = 0.03
    pipe.resample.output_shape = save_wcs['output_shape'].tolist()           #nx first and ny
    pipe.resample.crpix = save_wcs['crpix'].tolist()
    pipe.resample.crval = save_wcs['crval'].tolist()
    pipe.resample.rotation = save_wcs['rotation'].tolist()
    pipe.run(fits_files_corr)

    os.system('mv %sstep_i2d.fits %sstack_%s_pipe2_mywcs.fits'%(save_dir, save_dir, FILTER))
#################################################################################  
def main():
    #load reduction parameters
    if not os.path.isfile('reduction_params.yml'):
        sys.exit('no reduction_params.yml file found')

    with open('reduction_params.yml') as f:
    # use safe_load instead load
        dataMap = yaml.safe_load(f)

    os.environ["CRDS_CONTEXT"] = dataMap['pmap']
    basedir                    = dataMap['basedir']

    #run gaia align on
    gaia_align(basedir, 'F356W')

    #now stack this for reference
    restack(dict(FILTER='F356W', basedir=basedir))

    FILTERS = ['F200W','F115W']
    for FILTER in FILTERS:
        internal_align(basedir, FILTER)

    #now stack
    # with Pool(2) as pool:
    #     dictlist = [dict(FILTER=FILT, basedir=basedir) for FILT in FILTERS]
    #     pool.map(restack, dictlist)
##########################################################  
if __name__ == "__main__":
    main()
