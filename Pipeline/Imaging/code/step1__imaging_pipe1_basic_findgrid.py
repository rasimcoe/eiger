from glob import glob
import os
import shutil
import urllib
import sys
# Third Party Imports
from astropy.io import ascii as asc
from astropy.io import fits
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
import jwst
from jwst.datamodels import dqflags
from jwst.pipeline import calwebb_detector1
from jwst.pipeline import calwebb_image2
from jwst.pipeline import Image3Pipeline
from jwst import datamodels
from multiprocessing import Pool
from functools import partial
from jwst import resample
from jwst.resample import resample_utils
from jwst import datamodels
from astropy.modeling import models as astmodels
import yaml
from datetime import datetime
################################################
def detcal1(params):
    uncal_file = params['filename']
    workdir = params['workdir']

    detector1 = calwebb_detector1.Detector1Pipeline()
    detector1.output_dir = workdir
    detector1.save_results = True
    run_output = detector1.run(uncal_file)
    print(workdir, uncal_file) 
################################################
def detcal2(params):
    rate_file = params['filename']
    workdir = params['workdir']
    image2 = calwebb_image2.Image2Pipeline()
    image2.output_dir = workdir
    image2.save_results = True
    image2.resample.skip = True
    image2.run(rate_file)
################################################
def find_wcs_grid(all_cals, basedir):
    wcs_pickle = basedir + 'output_grid_wcs.npz'
    if not os.path.isfile(wcs_pickle):
        #load data models

        asn_exptypes = ['science']
        images = datamodels.open(all_cals, asn_exptypes=asn_exptypes)

        #get output wcs
        outwcs = resample_utils.make_output_wcs(images, pscale=0.03/3600, rotation=None)
        outhdr, junk = outwcs.to_fits()

        rotation = np.rad2deg(np.arctan2(outhdr['PC1_2'] , outhdr['PC2_2'])) % 360

        init_output_shape = np.array([outhdr['NAXIS1'], outhdr['NAXIS2']])
        init_crpix = np.array([outhdr['CRPIX1'], outhdr['CRPIX2']])
        crval = np.array([outhdr['CRVAL1'], outhdr['CRVAL2']])

        #slightly enlarge output frame
        dilate_output_shape = np.array([np.ceil(init_output_shape[0]/200)*200,
                               np.ceil(init_output_shape[1]/500)*500])
        dilate_output_shape = dilate_output_shape.astype(int)

        dilate_crpix = init_crpix + (dilate_output_shape-init_output_shape)/2.0


        #save 
        np.savez(wcs_pickle, rotation=rotation, output_shape=dilate_output_shape, crpix=dilate_crpix, crval=crval)
        return rotation, dilate_output_shape, dilate_crpix, crval
    else:
        save_wcs = np.load(wcs_pickle)
        return save_wcs['rotation'], save_wcs['output_shape'], save_wcs['crpix'], save_wcs['crval']
#################################################################################
def stack(params):
    FILTER = params['FILTER']
    basedir = params['basedir']

    wcs_pickle = basedir + 'output_grid_wcs.npz'
    save_wcs = np.load(wcs_pickle)

    workdir = basedir+FILTER +'/pipe1_basic/'
    fits_files_cal=glob(workdir + 'jw*_cal.fits')

    pipe = Image3Pipeline()
    pipe.output_dir=workdir
    pipe.tweakreg.skip=False
    pipe.tweakreg.align_to_gaia=True
    pipe.tweakreg.expand_refcat=True
    #pipe.tweakreg.min_gaia=4
    pipe.outlier_detection.skip=False
    pipe.save_results=True
    pipe.skymatch.skip=True
    #fixed grid for resample
    pipe.resample.pixel_scale = 0.03
    pipe.resample.output_shape = save_wcs['output_shape'].tolist()           #nx first and ny
    pipe.resample.crpix = save_wcs['crpix'].tolist()
    pipe.resample.crval = save_wcs['crval'].tolist()
    pipe.resample.rotation = save_wcs['rotation'].tolist()
    pipe.run(fits_files_cal)

    os.system('mv %sstep_i2d.fits %sstack_%s_pipe1.fits'%(workdir, workdir, FILTER))
################################################
def main():
    #load reduction parameters
    if not os.path.isfile('reduction_params.yml'):
        sys.exit('no reduction_params.yml file not found')

    with open('reduction_params.yml') as f:
    # use safe_load instead load
        dataMap = yaml.safe_load(f)

    os.environ["CRDS_SERVER_URL"]   = "https://jwst-crds.stsci.edu"
    os.environ["CRDS_DATA"]         = dataMap['cache']
    os.environ["CRDS_PATH"]         = dataMap['cache']
    os.environ["CRDS_CONTEXT"]      = dataMap['pmap']
    basedir                         = dataMap['basedir']


    n_procs = 40
    FILTERS = ['F115W', 'F200W', 'F356W']
    for FILTER in FILTERS:
        nseries = 8
        if FILTER == 'F356W': nseries = 2
        print(FILTER)
        download_dir = basedir + 'download/organised_output/IMAGING_'+FILTER+'/'
        workdir = basedir+FILTER +'/pipe1_basic/'
        os.chdir(workdir)

        #STEP 1
        fits_files_a=glob(download_dir + 'jw*_uncal.fits')
        #run the first 8 in series for crds to download the refence files, otherwise the parallel stalls
        for fname in fits_files_a[0:nseries]: 
            detcal1(dict(filename=fname, workdir=workdir))
        #now rest parallel
        with Pool(n_procs) as pool:
            dictlist = [dict(filename=fname, workdir=workdir) for fname in fits_files_a[nseries:]]
            pool.map(detcal1, dictlist)

        #STEP 2
        fits_files_rate=glob(workdir + 'jw*_rate.fits')
        #run the first 8 in series for crds
        for fname in fits_files_rate[0:nseries]: 
            detcal2(dict(filename=fname, workdir=workdir))
        #now rest parallel
        with Pool(n_procs) as pool:
            dictlist = [dict(filename=fname, workdir=workdir) for fname in fits_files_rate[nseries:]]
            pool.map(detcal2, dictlist)

    #find missing cal files
    for FILTER in FILTERS:
        missed_list = []
        missed_rate = []
        download_dir = basedir + 'download/organised_output/IMAGING_'+FILTER+'/'
        workdir = basedir+FILTER +'/pipe1_basic/'
        all_uncal = glob(download_dir + 'jw*uncal.fits')
        for fname in all_uncal:
            calname = workdir + (fname.split('/')[-1]).replace('uncal.','cal.')
            if not os.path.isfile(calname):
                print('missing file ', FILTER, fname.split('/')[-1])
                missed_list.append(fname)
                missed_rate.append(workdir + (fname.split('/')[-1]).replace('uncal.','rate.'))

        if len(missed_list) > 0:
            with Pool(n_procs) as pool:
                dictlist = [dict(filename=fname, workdir=workdir) for fname in missed_list]
                pool.map(detcal1, dictlist)
            #now rest parallel
            with Pool(n_procs) as pool:
                dictlist = [dict(filename=fname, workdir=workdir) for fname in missed_rate]
                pool.map(detcal2, dictlist)     
    
    #work out grid
    all_cals = glob(basedir + '*/pipe1_basic/jw*_cal.fits')
    find_wcs_grid(all_cals, basedir)
    
    #now stack in parallel
    # with Pool(3) as pool:
    #     dictlist = [dict(FILTER=FILT, basedir=basedir) for FILT in FILTERS]
    #     pool.map(stack, dictlist)

    stack(dict(FILTER='F356W', basedir=basedir))
    stack(dict(FILTER='F200W', basedir=basedir))
    stack(dict(FILTER='F115W', basedir=basedir))
##########################################################  
if __name__ == "__main__":
    main()
