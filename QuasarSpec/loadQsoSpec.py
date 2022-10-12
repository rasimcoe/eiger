from eiger.Database.SQL import eigerdb
from eiger.Database.Fileserver import s3_etag
import boto3
from astropy.io import fits
import os
import numpy as np

#############################################################################
#
# Support routine to check local and remote hashes of a file in S3
# and returns either the local file (as appropriate) or downloads
# the remote one to the cache and returns the result.
#

def getSpec(specfile, offline=False, filenames=False):

    localfile = os.getenv('EIGER_CACHE')+'/'+specfile[1]

    if (not offline):
        try:

            # The etag is a MD5 hash stored by AWS for each unique file
            remote_etag = s3_etag.etag_remotehash(specfile[0],specfile[1])

            if (os.path.exists(localfile)):        
                local_etag  = s3_etag.etag_localhash(localfile)
            else:
                local_etag  = ''
                
            if (remote_etag != local_etag):

                # If the remote and local hash are out of sync, it means
                # there is a newer version uploaded to the cloud.  Go and
                # get it.
                
                # Make sure that the local directories exist to write to Cache
                os.makedirs(os.path.dirname(localfile),exist_ok=True)
                
                # Get the appropriate S3 credentials
                s3_resource = boto3.resource('s3')
                reply = s3_resource.Bucket(specfile[0]).download_file(\
                                                    specfile[1],localfile)
                print("   Downloading remote version of spectrum from the AWS cloud")

            else:
                # If the hashes align, then the local cache has
                # an up-to-date version 
                print("   Using local cached spectrum")

        except:
            print("ERROR: AWS file could not be retrieved")
            
    return(localfile)


#############################################################################

def parseSpec(fitsfile, instrument):

    outspec = {}

    np.seterr(divide='ignore', invalid='ignore')
    
    if (instrument == 'FIRE'):
        tmp = fits.open(fitsfile)[1].data
        if ('ivar' in tmp.columns.names):
            # mask = tmp['mask']
            mask = np.array(tmp['mask'],dtype=bool)
            outspec['wave'] = tmp['wave'][mask]
            outspec['flux'] = tmp['flux'][mask]
            outspec['ivar'] = tmp['ivar'][mask]
        else:
            outspec['wave'] = tmp['wave'][0]
            outspec['flux'] = tmp['flux'][0]
            outspec['ivar'] = 1.0/tmp['sig'][0]**2
            outspec['ivar'][np.isinf(outspec['ivar'])] = 0.0
            
    elif (instrument == 'HIRES'):
        tmp = fits.open(fitsfile)[1].data[0]
        outspec['wave'] = tmp['wave']
        outspec['flux'] = tmp['flux']
        outspec['ivar'] = 1.0/tmp['sig']**2
        outspec['ivar'][np.isinf(outspec['ivar'])] = 0.0
        outspec['mask'] = np.ones(len(tmp['wave']))
        
    elif (instrument == 'XShooter'):
        tmp = fits.open(fitsfile)[1].data
        mask = np.array(tmp['mask'],dtype=bool)
        outspec['wave'] = tmp['wave'][mask]
        outspec['flux'] = tmp['flux'][mask]
        outspec['ivar'] = tmp['ivar'][mask]
        
    elif (instrument == 'MOSFIRE'):
        tmp = fits.open(fitsfile)[1].data
        mask = np.array(tmp['mask'],dtype=bool)
        # mask = tmp['mask']
        outspec['wave'] = tmp['wave'][mask]
        outspec['flux'] = tmp['flux'][mask]
        outspec['ivar'] = tmp['ivar'][mask]

    elif (instrument == 'FIRE_XSH'):
        tmp = fits.open(fitsfile)[1].data
        mask = np.array(tmp['mask'],dtype=bool)
        outspec['wave'] = tmp['wave'][mask]
        outspec['flux'] = tmp['flux'][mask]
        outspec['ivar'] = tmp['ivar'][mask]
        
    contname = fitsfile[:-5]+'_contin.fits'
    if(os.path.exists(contname)):
        tmp = fits.open(contname)[1].data
        try:
            outspec['cont'] = tmp['cont'][mask]
        except:
            outspec['cont'] = tmp['cont']
    else:
        # Go and get it from S3
        try:
            cc = contname.split('//')[1]
            print(f"Fetching continuum file from S3 cloud ({cc})")
            getSpec(['gto1243',cc])
            tmp = fits.open(contname)[1].data
            try:
                outspec['cont'] = tmp['cont'][mask]
            except:
                outspec['cont'] = tmp['cont']
            # outspec['inlier_mask'] = tmp['mask']
        except:
            print("WARNING: No valid continuum spectrum exists for this object, either locally or on the cloud.")
            
    np.seterr(divide='warn', invalid='warn')

    return(outspec)
        
#############################################################################
#############################################################################
#############################################################################
#############################################################################



def loadQsoSpec(obj_id, spectrographs=['XShooter', 'FIRE', 'MOSFIRE', 'HIRES','FIRE_XSH'], \
                revision='current', offline=False):

    if (obj_id < 1 or obj_id > 6):
        print("ERROR: Quasar ID must be between 1 and 6")
        return(None)
    
    spectra = {}

    if (offline == False):
        ##### Find the spectra that exist in the observations database
        db = eigerdb.Eigerdb()
        db.getcursor()

        reply = db.query(f"select name,zem from Quasars where id={obj_id}")
        print(f"Object: {reply[0][0]}, quasarid={obj_id}")
        spectra['objid']   = obj_id
        spectra['objname'] = reply[0][0]
        spectra['z_em'] = reply[0][1]
        spectra['revision'] = revision
    else:
        if (obj_id == 4):
            spectra['objid']   = obj_id
            spectra['objname'] = 'J0100'
            spectra['z_em'] = 6.3258
            spectra['revision'] = 'current'
        else:
            print(f"loadQsoSpec not configured for objid == {obj_id}")
            return()
    
    for spectrograph in spectrographs:

        if (offline == False):
            querystring = """
            SELECT awsbucket,awspath 
            FROM Observations 
            WHERE quasarid={} 
            AND revision=\'{}\'
            AND instrument=\'{}\'""".format(obj_id,revision,spectrograph)
            
            obs = db.query(querystring)

        else:
            if (obj_id == 4):
                if (spectrograph=='XShooter'):
                    obs = [['gto1243','VLT/XShooter/J0100+2802_VIS_XSHOOTER.fits'], \
                           ['gto1243','VLT/XShooter/J0100+2802_NIR_XSHOOTER.fits']]
                elif (spectrograph=='FIRE'):
                    obs = [['gto1243','Magellan/FIRE/J0100+28_FIRE_firehose.fits']]
                elif (spectrograph == 'HIRES'):
                    obs = [['gto1243','Keck/HIRES/J0100_HIRES.fits']]
                else:
                    print("ERROR: bad instrument passed in offline mode")
                    return()
                    
        if (len(obs) > 0):

            print(f"{spectrograph}:")

            if (spectrograph == 'XShooter'):
                # Grab both the VIS and NIR arms
                for oo in obs:
                    local_file = getSpec(oo, offline=offline)
                    if ('VIS' in local_file):
                        arm = 'XSH_VIS'
                    else:
                        arm = 'XSH_NIR'

                    spectrum = parseSpec(local_file,spectrograph)
                    spectra[arm] = spectrum

            elif (spectrograph == 'MOSFIRE'):
                # Separate files for Y, J, H, K
                for oo in obs:
                    local_file = getSpec(oo, offline=offline)
                    if ('_Y_' in local_file):
                        arm='MOSFIRE_Y'
                    elif ('_J_' in local_file):
                        arm='MOSFIRE_J'
                    elif ('_H_' in local_file):
                        arm='MOSFIRE_H'
                    elif ('_K_' in local_file):
                        arm='MOSFIRE_K'

                    spectrum = parseSpec(local_file,spectrograph)
                    spectra[arm] = spectrum

            else:
                local_file = getSpec(obs[0], offline=offline)
                spectrum = parseSpec(local_file,spectrograph)
                spectra[spectrograph] = spectrum

        else:
            continue
            
    if (offline == False):
        db.close()
    
    return(spectra)

