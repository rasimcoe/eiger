import os
from glob import glob
import shutil
from astropy.io import fits
from astroquery.mast import Observations
import sys
import argparse
import getpass
import yaml
##########################################################
def check_complete(savedir, science_list):
    #check complete
    #science_list = my_session.query_criteria(project='JWST', proposal_id='1243', calib_level=2)
    for sci in science_list:
        newname = savedir +'*/' + sci['obs_id'] + '_uncal.fits'
        filelist = glob(newname)
        if len(filelist) == 0:
            print('missing uncal file ',sci['obs_id'])
##########################################################
def check_EIGER_complete(savedir):
    filelist = glob(savedir +'IMAGING_F115W/jw01243*_uncal.fits')
    obs_base = filelist[0].split('/')[-1][0:10]
    for fname in filelist:
        if fname.split('/')[-1][0:10] != obs_base: print('more than one obs id!')

    #loop round main observation structure
    obs_str = [ dict(FILTER='F115W', seqid='02101', nint=12, nvis=4, modules=['a','b'], cams=['1','2','3','4']), #grism a
                dict(FILTER='F200W', seqid='02103', nint=12, nvis=4, modules=['a','b'], cams=['1','2','3','4']), #grism b
                dict(FILTER='F200W', seqid='02104', nint=1 , nvis=4, modules=['a','b'], cams=['1','2','3','4']), #direct
                dict(FILTER='F200W', seqid='02106', nint=2 , nvis=4, modules=['a','b'], cams=['1','2','3','4']), #out of field
                dict(FILTER='F356W', seqid='02104', nint=1 , nvis=4, modules=['a','b'], cams=['long']),          #direct
                dict(FILTER='F356W', seqid='02106', nint=2 , nvis=4, modules=['a','b'], cams=['long'])]          #out of field    
    for group in obs_str:
        for iint in range(1,group['nint']+1):
            for ivis in range(1,group['nvis']+1):
                for mod in group['modules']:
                    for cam in group['cams']:
                        filename = '%s%03i_%s_%05i_nrc%s%s_uncal.fits' % (obs_base, ivis, group['seqid'], iint, mod, cam)
                        if not os.path.isfile(savedir +'IMAGING_'+group['FILTER']+'/'+filename):
                            print('%s missing file, main obs %s' % (filename, group['FILTER']))

    filelist = glob(savedir +'IMAGING_F356W/jw01243*02101*_uncal.fits')
    ext_base = filelist[0].split('/')[-1][0:10]
    for fname in filelist:
        if fname.split('/')[-1][0:10] != ext_base: print('more than one obs id!')

    ext_std = [dict(FILTER='F200W', seqid=['02101','02101'], nint=1 , nvis=2, modules=['a','b'], cams=['1','2','3','4']), #extra
               dict(FILTER='F356W', seqid=['02101','02101'], nint=1 , nvis=2, modules=['a','b'], cams=['long'])]          #extra                      
    
    for group in ext_std:
        for iint in range(1,group['nint']+1):
            for ivis in range(1,group['nvis']+1):
                for mod in group['modules']:
                    for cam in group['cams']:
                        filename = '%s%03i_%s_%05i_nrc%s%s_uncal.fits' % (ext_base, ivis, group['seqid'][ivis-1], iint, mod, cam)
                        if not os.path.isfile(savedir +'IMAGING_'+group['FILTER']+'/'+filename):
                            print('%s missing file, extra obs %s' % (filename, group['FILTER']))

    # token = getpass.getpass(prompt='MAST token: ', stream=None)
    # my_session = Observations(mast_token=token)
    # science_list = my_session.query_criteria(project='JWST', proposal_id='1243', target_name=qso_name, dataproduct_type='image')
    # tab = my_session.download_products(test, mrp_only=False, download_dir='.',  productSubGroupDescription=['UNCAL'])
    # tab['Message']
##########################################################
def organise_output(workdir, savedir):
    file_list = glob(workdir + 'mastDownload/JWST/jw012430*/*uncal.fits')

    for filename in file_list:
        hdr = fits.getheader(filename)
        if hdr['EXP_TYPE'] == 'NRC_WFSS':
            mode = 'GRISM_'
        else:
            mode = 'IMAGING_'

        modedir = savedir + mode + hdr['FILTER'] + '/'
        newname = modedir + filename.split('/')[-1]
        if not os.path.isfile(newname): 
            if not os.path.isdir(modedir): os.mkdir(modedir)
            shutil.copy(filename, newname)
##########################################################
def MAST_download(qso_name):
    #open session
    token = getpass.getpass(prompt='MAST token: ', stream=None)
    my_session = Observations(mast_token=token)
    #sessioninfo = my_session.session_info()
    #fetch filelist
    science_list = my_session.query_criteria(project='JWST', proposal_id='1243', target_name=qso_name, dataproduct_type='image')
    print(science_list)
    if len(science_list) == 0:
        sys.exit("MAST query returned no rows, check target name and token")
    #download uncal files
    tab = my_session.download_products(science_list['obsid'], mrp_only=False, download_dir='.',  productSubGroupDescription=['UNCAL'])
    tab['Message']
    return science_list
##########################################################  
def make_directories(basedir):
    FILTERS = ['F115W', 'F200W', 'F356W']
    for FILTER in FILTERS:
        filtdir = basedir+FILTER+'/'
        if not os.path.isdir(filtdir): os.mkdir(filtdir)

        for dir_type in ['pipe1_basic', 'pipe2_mywcs', 'pipe3_skyfix', 'mycals', 'pipe4_filt']:
            workdir = filtdir+dir_type+'/'
            if not os.path.isdir(workdir): os.mkdir(workdir)
##########################################################  
def write_paramfile(qso_name, basedir, pmap, cache):
    data =  {'qso_name': qso_name,
                         'basedir': basedir,
                         'pmap': pmap,
                         'cache': cache}

    with open('reduction_params.yml', 'w') as yaml_file:
        yaml.dump(data, yaml_file, default_flow_style=False)

    # with open('reduction_params.yml') as f:
    # # use safe_load instead load
    #     dataMap = yaml.safe_load(f)
##########################################################    

def main():
    p = argparse.ArgumentParser()
    p.add_argument("-q", "--qso", help='quasar name on MAST  e.g. "2MASS J01001301+2802257" "QSO J1120+0641" ')
    p.add_argument("-d", "--dir", help='base directory e.g. "/scratch/mruari/EIGER/imaging/J0100+2802/"')
    p.add_argument("-p", "--pmap", help='CRDS pmap e.g. "jwst_0988.pmap"')
    p.add_argument("-c", "--cache", help='CRDS cache directory e.g. "/scratch/mruari/EIGER/cache/crds_cache"')

    args = p.parse_args() 

    #QSO_name = 'J0100+2802'
    #basedir = '/scratch/mruari/EIGER/imaging/'+QSO_name+'/'

    workdir = args.dir + 'download/'
    savedir = args.dir + 'download/organised_output/'

    print( args.dir, args.qso )
    if not os.path.isfile('reduction_params.yml'):
        write_paramfile(args.qso, args.dir, args.pmap, args.cache)
    #if os.path.isfile('reduction_params.yml'):
    #    sys.exit('existing reduction_params.yml found: delete or rename')

    if not os.path.isdir(workdir): os.mkdir(workdir)
    os.chdir(workdir)
    science_list = MAST_download(args.qso)


    if not os.path.isdir(savedir): os.mkdir(savedir)
    organise_output(workdir, savedir)

    make_directories(args.dir)

    #check_complete(savedir, science_list)
    check_EIGER_complete(savedir)
##########################################################  
if __name__ == "__main__":
    main()