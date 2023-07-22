from astroquery.mast import Observations
#from ccdproc import ImageFileCollection
from astropy.table import Table
from astropy.io import fits
from astropy.io import ascii
from astropy.visualization import ZScaleInterval
from IPython.display import Image
from glob import glob
import matplotlib.pyplot as plt
import numpy as np
import os
import shutil
from drizzlepac import tweakreg
from drizzlepac import astrodrizzle
from astropy.nddata import NDData
from photutils import detect_threshold, DAOStarFinder, IRAFStarFinder
from astropy.table import Table, Column
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
from photutils.psf import extract_stars
from astroquery.gaia import Gaia
import yaml
from astropy.time import Time
import sys
from stwcs.wcsutil import HSTWCS
from tweakwcs import fit_wcs, align_wcs, FITSWCSCorrector, XYXYMatch, FITSWCS
from drizzlepac import updatehdr
from acstools import acszpt
from astropy.time import Time
#####################################################################################
def sexcat(filename, name, wcs):
    paramfile = 'config.ast'

    os.system('cp /net/galaxy-data/export/galaxydata/mruari/EIGER/imaging/param_files/astrometry_check_config/* .')

    catname_det = 'astcat_'+name+'.fits'
    oargs = '-CATALOG_NAME %s' % (catname_det) 
    cmd = 'source-extractor '+filename+'  -c '+paramfile+' '+oargs
    os.system(cmd)

    cat = Table.read(catname_det)
    cat.rename_column('ALPHA_J2000','RA')
    cat.rename_column('DELTA_J2000','DEC')
    cat.rename_column('X_IMAGE','x')
    cat.rename_column('Y_IMAGE','y')

    cat['RA'], cat['DEC'] = wcs.all_pix2world(cat['x'], cat['y'],0)
    return cat
#####################################################################################
def tweak_nrcref(refimg, blindstack, ir=False):

    if ir:
        exten = 'flt.fits'
        imagefindcfg = {'threshold': 5, 'conv_width': 2.0}
    else:
        exten = 'flc.fits'
        imagefindcfg = {'threshold': 500, 'conv_width': 3.5, 'dqbits': ~4096}

    #load images
    temp_stack = fits.open(blindstack)
    ref_stack = fits.open(refimg)

    #open wcs
    refwcs = HSTWCS(ref_stack, ('SCI', 1))
    tempwcs = HSTWCS(temp_stack)

    #make catalogs
    refsex = sexcat(refimg+'[1]',    "F115W_ref",  refwcs)
    imsex  = sexcat(blindstack+'[0]', "blind_temp", refwcs)
    imwcs = FITSWCSCorrector(tempwcs)

    # Match sources in the catalogs:
    match = XYXYMatch(searchrad=10, separation=0.01, tolerance=5, use2dhist=True)
    ridx, iidx = match(refsex, imsex, imwcs)

    # Align image WCS:
    fit_res = fit_wcs(refsex[ridx], imsex[iidx], imwcs,fitgeom='rshift')
    aligned_imwcs = fit_res.wcs
    #imcat.meta['aligned_wcs'] = aligned_imwcs

    #propigate tweak to flc.fits
    reftwc = FITSWCS(refwcs)
    flcs = glob('hst*_flc.fits')
    for fname in flcs:
        hdul = fits.open(fname)

        flcwcs = HSTWCS(hdul, ('SCI', 1))
        tpwcs_corrector = FITSWCSCorrector(flcwcs)    
        tpwcs_corrector.set_correction(matrix=fit_res.meta['matrix'], shift=-fit_res.meta['shift'], ref_tpwcs=reftwc)
        #tpwcs_corrector.set_correction(shift=[-10,0], ref_tpwcs=reftwc)
        #update header
        updatehdr.update_wcs(hdul, 1, tpwcs_corrector.wcs, wcsname='TWEAK', reusename=True, verbose=False)

        flcwcs2 = HSTWCS(hdul, ('SCI', 2))
        tpwcs_corrector = FITSWCSCorrector(flcwcs2)    
        tpwcs_corrector.set_correction(matrix=fit_res.meta['matrix'], shift=-fit_res.meta['shift'], ref_tpwcs=reftwc)
        #tpwcs_corrector.set_correction(shift=[-10,0], ref_tpwcs=reftwc)
        #update header
        updatehdr.update_wcs(hdul, 4, tpwcs_corrector.wcs, wcsname='TWEAK', reusename=True, verbose=False)

        hdul.writeto(fname, overwrite=True)

    #restack
    stackname = FILTER+'_stack_gaiaref_nrctweak_final_skylocalmin2'

    astrodrizzle.AstroDrizzle('hst_*%s'%(exten),
              output=stackname,
              preserve=False,
              clean=False,
              build=False,
              context=False,
              driz_sep_bits='256,64,16',
              combine_type='minmed',
              final_bits='256,64,16',
              runfile='drizzle.log',
              skymethod='localmin', 
              driz_separate=True,
              median=True,
              blot=True,
              driz_cr=True,
              final_wcs=True,
              final_refimage=refimg)
#####################################################################################
def get_zpts(stackname, FILTER):
    hdu = fits.open(stackname)

    mean_mjd = (hdu[0].header['EXPSTART']+hdu[0].header['EXPEND'])/2. 
    exp_time = Time(mean_mjd, format='mjd')

    q = acszpt.Query(date=exp_time.isot[:10], detector="WFC", filt=FILTER)
    zpt_table = q.fetch()

    fits.setval(stackname, 'PHOTFLAM', value=zpt_table['PHOTPLAM'].value[0], ext=0)
    fits.setval(stackname, 'ZPT_AB', value=zpt_table['ABmag'][0].value, ext=0)
#####################################################################################
def obj_q( qso_name):
    # Retrieve the IACS FLC data prodcuts
    #MAST_names are target names in relevent data 
    obs_table = Observations.query_object(qso_name)

    obs_sel = obs_table[(obs_table['project'] == 'HST')&(obs_table['distance'] < 300)]
    opt = np.unique(obs_sel[(obs_sel['instrument_name'] == 'WFC3/UVIS')|(obs_sel['instrument_name'] == 'ACS/WFC')]['filters']).data.data
    ir  = np.unique(obs_sel[(obs_sel['instrument_name'] == 'WFC3/IR')]['filters']).data.data
    #filters =  np.unique(obs_sel['filters']).data
    mast_names = np.unique( obs_sel['target_name']).data.data
    for io in range(len(opt)):
        opt[io] = (opt[io].replace('CLEAR1L;','')).replace(';CLEAR2L','')
    return opt, ir, mast_names
#####################################################################################
def fetch_data(FILTER, mast_names, ir=False):
    # Retrieve the IACS FLC data prodcutsobs_sel['instrument']
    #mast_names are target names in relevent data 
    science_list = Observations.query_criteria(filters=FILTER, project='HST', target_name=mast_names)

    if ir:
        Observations.download_products(science_list['obsid'], mrp_only=False, download_dir='.', productSubGroupDescription=['FLT'])
    else:
        Observations.download_products(science_list['obsid'], mrp_only=False, download_dir='.', productSubGroupDescription=['FLC'])
    science_files = glob(os.path.join(os.curdir, 'mastDownload', 'HST', 'hst*', '*fits'))
    for filename in science_files:
        shutil.copy(filename, filename.split('/')[-1])

    # collect_img = ImageFileCollection('./', glob_include="*flc.fits", ext=0,
    #                                  keywords=["asn_id", "detector", "filter", "nsamp",
    #                                            "exptime", "postarg1", "postarg2"])

    # exp_table = collect_img.summary
    # exp_table['exptime'].format = '7.1f'
    # exp_table['postarg1'].format = '7.2f'
    # exp_table['postarg2'].format = '7.2f'
    # exp_table
#####################################################################################
def drizzle_flcs(FILTER, ir=False):
    if ir:
        astrodrizzle.AstroDrizzle('*flt.fits',
                      output=FILTER+'_stack_F356Wref_final',
                      preserve=False,
                      clean=True,
                      build=False,
                      context=False,
                      driz_sep_bits='256,64,16',
                      combine_type='minmed',
                      final_bits='256,64,16',
                      runfile=FILTER+'_blind.log',
                      skymethod='localmin', 
                      driz_separate=True,
                      median=True,
                      blot=True,
                      driz_cr=True,
                      final_wcs=True,
                      final_refimage=basedir + 'F356W/pipe2_mywcs/stack_F356W_pipe2_mywcs.fits')
    else:
        astrodrizzle.AstroDrizzle('*flc.fits',
                      output=FILTER+'_stack_F356Wref_final_skylocalmin',
                      preserve=False,
                      clean=True,
                      build=False,
                      context=False,
                      driz_sep_bits='256,64,16',
                      combine_type='minmed',
                      final_bits='256,64,16',
                      runfile=FILTER+'_blind.log',
                      skymethod='localmin', 
                      driz_separate=True,
                      median=True,
                      blot=True,
                      driz_cr=True,
                      final_wcs=True,
                      final_refimage=basedir + 'F356W/pipe2_mywcs/stack_F356W_pipe2_mywcs.fits')
#####################################################################################
def make_star_cat(SCI, WHT, header, FILTER):
    
    mask = (WHT == 0)|(SCI == 0)

    mean, median, std = sigma_clipped_stats(SCI, sigma=3.0, mask=mask)  

    daofind = DAOStarFinder(fwhm=1.8, threshold=200.*std,
                            sharplo=0.5, sharphi=1.0, roundlo=-1.0,
                            roundhi=1.0, brightest=None,
                            peakmax=100, exclude_border=True)

    sources = daofind(SCI - median, mask=mask)
    sources.rename_column('xcentroid','x')
    sources.rename_column('ycentroid','y')

    #convert to ra dec
    wcs = WCS(header)
    ra, dec = wcs.all_pix2world(sources['x'], sources['y'],0)
    sources.add_column(Column(ra, name='RA'))
    sources.add_column(Column(dec, name='DEC'))

    if FILTER == 'F356W':
        nddata = NDData(SCI)  
        all_stars = extract_stars(nddata, sources, size=31) 
        saturated = np.any(np.array(all_stars.data) == 0, axis=(1,2))
        clean_sources = sources[(~saturated)]
    else:
        #remove sources near edges
        hsize = 20
        mask = ((sources['x'] > hsize) & (sources['x'] < (SCI.shape[1] -1 - hsize)) &
         (sources['y'] > hsize) & (sources['y'] < (SCI.shape[0] -1 - hsize)))  
        sel_sources = sources[(mask)]
        nddata = NDData(WHT)  
        all_stars = extract_stars(nddata, sel_sources, size=31) 
        saturated = np.any(np.array(all_stars.data) == 0, axis=(1,2))
        clean_sources = sel_sources[(~saturated)]
    #remove close dups
    #keep_sources = isolated(sources, 2.0)
    columns = ['RA', 'DEC', 'id', 'flux']
    catalog = clean_sources[columns]
    clean_sources.write('star_cat_'+FILTER+'.fits', overwrite=True)
    #sources.write('ref_cat_F356W_v4.txt', format='ascii.commented_header', overwrite=True)
    return clean_sources
###############################################################################
def internal_tweakreg(FILTER, filtdir, basedir, ir=False):

    ref_filename = basedir + 'F356W/pipe2_mywcs/stack_F356W_pipe2_mywcs.fits'
    ref_img = fits.open(ref_filename)
    #make align cat
    ref_cat = make_star_cat(ref_img['SCI'].data, ref_img['WHT'].data, ref_img['SCI'].header, 'F356W')
    #write refcat
    ref_cat.rename_column('RA', 'ra')
    ref_cat.rename_column('DEC', 'dec')
    reduced_query = ref_cat['ra', 'dec', 'flux']
    reduced_query.write('F356W.cat', format='ascii.commented_header', overwrite=True)

    refcat = 'F356W.cat'
    wcsname ='F356W_ref'

    if ir:
        tweakreg.TweakReg('*flt.fits',  # Pass input images
            updatehdr=True,  # update header with new WCS solution
            imagefindcfg={'threshold': 5, 'conv_width': 2.0},  # Detection parameters, threshold varies for different data
            refcat=refcat,  # Use user supplied catalog (Gaia)
            interactive=False,
            see2dplot=False,
            shiftfile=True,  # Save out shift file (so we can look at shifts later)
            outshifts='F356W_rscale.txt',  # name of the shift file
            wcsname=wcsname,  # Give our WCS a new name
            reusename=True,
            sigma=3.0,
            nclip=2,
            fitgeometry='rscale', # Use the 6 parameter fit
            minobj=3, 
            searchrad=10.0, 
            tolerance=5.,
            clean=True) 
    else:
        tweakreg.TweakReg('*flc.fits',  # Pass input images
            updatehdr=True,  # update header with new WCS solution
            imagefindcfg={'threshold': 1000, 'conv_width': 6.0, 'dqbits': ~4096},  # Detection parameters, threshold varies for different data
            refcat=None,  # Use user supplied catalog (Gaia)
            interactive=False,
            see2dplot=False,
            shiftfile=True,  # Save out shift file (so we can look at shifts later)
            outshifts='blind_rscale.txt',  # name of the shift file
            wcsname=refcat,  # Give our WCS a new name
            reusename=True,
            sigma=3.0,
            nclip=2,
            fitgeometry='rscale', # Use the 6 parameter fit
            minobj=3, 
            searchrad=10.0, 
            tolerance=5.,
            clean=True) 
#################################################################################
def Gaia_drizzle(FILTER, refimg, ir= False):

    if ir:
        exten = 'flt.fits'
        imagefindcfg = {'threshold': 5, 'conv_width': 2.0}
    else:
        exten = 'flc.fits'
        imagefindcfg = {'threshold': 500, 'conv_width': 3.5, 'dqbits': ~4096}

    all_flcs = glob('hst_*'+exten)
    test = fits.getheader(all_flcs[0])

    #now getch catalog
    g_rad = 10.0
    coord = SkyCoord(test['RA_TARG'], test['DEC_TARG'], unit=(u.deg, u.deg), frame='icrs')
    radius = u.Quantity(g_rad, u.arcmin)

    gaia_query = Gaia.query_object_async(coordinate=coord, radius=radius)

    #break into epocs
    prog_ids = np.unique([cur.split('_')[1] for cur in all_flcs])

    for prog in prog_ids:
        current_files = glob('hst_%s*%s'%(prog, exten))

        #get epoch
        example = fits.getheader(current_files[0])
        mjd_date = example['EXPSTART']
        epoch_obs = Time(mjd_date, format='mjd')
        #apply proper motions

        coords = SkyCoord(ra=gaia_query['ra'], dec=gaia_query['dec'], 
                         pm_ra_cosdec=gaia_query['pmra'],
                         pm_dec=gaia_query['pmdec'],
                         obstime=Time(gaia_query['ref_epoch'], format='jyear', scale='tcb'))

        
        c_obs_epoch = coords.apply_space_motion(epoch_obs)

        epoch_cat = Table(gaia_query['ra', 'dec', 'phot_g_mean_mag'])
        epoch_cat['ra'] = c_obs_epoch.ra.value
        epoch_cat['dec'] = c_obs_epoch.dec.value
        epoch_keep = epoch_cat[(~np.isnan(epoch_cat['ra']))&(~np.isnan(epoch_cat['dec']))]

        #save reduced catalog
        refcat = 'gaia_%s.cat' % (prog)
        cw = 3.5  # Set to two times the FWHM of the PSF.
        wcsname = 'Gaia'  # Specify the WCS name for this alignment
        thresh = 500.0  
        epoch_keep.write(refcat, format='ascii.commented_header')

        #apply tweaks
        tweakreg.TweakReg('hst_%s*%s'%(prog, exten),  # Pass input images
                        updatehdr=True,  # update header with new WCS solution
                        imagefindcfg=imagefindcfg,  # Detection parameters, threshold varies for different data
                        refcat=refcat,  # Use user supplied catalog (Gaia)
                        interactive=False,
                        see2dplot=False,
                        shiftfile=True,  # Save out shift file (so we can look at shifts later)
                        outshifts='Gaia_rscale_%s.txt'%(prog),  # name of the shift file
                        wcsname=wcsname,  # Give our WCS a new name
                        reusename=True,
                        sigma=3.0,
                        nclip=2,
                        fitgeometry='rscale', # Use the 6 parameter fit
                        minobj=1, 
                        searchrad=10.0, 
                        tolerance=5.)  
                        
        astrodrizzle.AstroDrizzle('hst_%s*%s'%(prog, exten),
                      output='stack_hst_prog_%s'%(prog),
                      preserve=False,
                      clean=False,
                      build=False,
                      context=False,
                      driz_sep_bits='256,64,16',
                      combine_type='minmed',
                      final_bits='256,64,16',
                      runfile=FILTER+'_blind.log',
                      skymethod='localmin', 
                      driz_separate=True,
                      median=True,
                      blot=True,
                      driz_cr=True,
                      final_wcs=True,
                  final_refimage=refimg)

                    
    astrodrizzle.AstroDrizzle('hst_*%s'%(exten),
                  output=FILTER+'_stack_gaiaref_final_skylocalmin2',
                  preserve=False,
                  clean=False,
                  build=False,
                  context=False,
                  driz_sep_bits='256,64,16',
                  combine_type='minmed',
                  final_bits='256,64,16',
                  runfile=FILTER+'_blind.log',
                  skymethod='localmin', 
                  driz_separate=True,
                  median=False,
                  blot=False,
                  driz_cr=False,
                  final_wcs=True,
                  final_refimage=refimg)
    return FILTER+'_stack_gaiaref_final_skylocalmin2'+'_drc_sci.fits'
#################################################################################
def main():
    #load reduction parameters
    if not os.path.isfile('reduction_params.yml'):
        sys.exit('no reduction_params.yml file found')

    with open('reduction_params.yml') as f:
    # use safe_load instead load
        dataMap = yaml.safe_load(f)

    qso_name                   = dataMap['qso_name']
    basedir                    = dataMap['basedir']

    refimg = asedir + 'F356W/pipe2_mywcs/stack_F356W_pipe2_mywcs.fits'


    if not os.path.isdir(basedir+'HST/'): os.mkdir(basedir+'HST/')


    #first basic query 
    opt_filters, ir_filters, mast_names = obj_q( qso_name)

    for FILTER in opt_filters:
        filtdir = basedir+'HST/'+FILTER+'/'
        if not os.path.isdir(filtdir): os.mkdir(filtdir)
        os.chdir(filtdir)

        #download data
        fetch_data(FILTER, mast_names)

        #tweak reg to Gaia
        step1 = Gaia_drizzle(FILTER, refimg)

        #nrc tweak and restack
        tweak_nrcref(refimg, step1)

        #add ZPT 
        get_zpts(FILTER+'_stack_gaiaref_nrctweak_final_skylocalmin2'+'_drc_sci.fits', FILTER)

    for FILTER in ir_filters:
        filtdir = basedir+'HST/'+FILTER+'/'
        if not os.path.isdir(filtdir): os.mkdir(filtdir)
        os.chdir(filtdir)

        #download data
        fetch_data(FILTER, mast_names, ir=True)

        #tweak reg to Gaia
        step1 = Gaia_drizzle(FILTER, refimg, ir=True)

        #nrc tweak and restack
        tweak_nrcref(refimg, step1, ir=True)

##########################################################  
if __name__ == "__main__":
    main()