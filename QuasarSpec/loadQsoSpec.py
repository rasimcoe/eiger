from eiger.Database.SQL import eigerdb
from eiger.Database.Fileserver import s3_etag
import boto3
from astropy.io import fits
import os
from numpy import sqrt

#############################################################################
#
# Support routine to check local and remote hashes of a file in S3
# and returns either the local file (as appropriate) or downloads
# the remote one to the cache and returns the result.
#

def getSpec(specfile):
    
    try:

        # The etag is a MD5 hash stored by AWS for each unique file
        remote_etag = s3_etag.etag_remotehash(specfile[0],specfile[1])

        localfile = os.getenv('EIGER_CACHE')+'/'+specfile[1]
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
            reply = s3_resource.Bucket(specfile[0]).download_file(specfile[1],localfile)
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
    
    if (instrument == 'FIRE'):
        tmp = fits.open(fitsfile)[1].data[0]
        outspec['wave'] = tmp['wave']
        outspec['flux'] = tmp['flux']
        outspec['sig'] = tmp['sig']

    elif (instrument == 'HIRES'):
        tmp = fits.open(fitsfile)[1].data[0]
        outspec['wave'] = tmp['wave']
        outspec['flux'] = tmp['flux']
        outspec['sig'] = tmp['sig']
        
    elif (instrument == 'XShooter'):
        tmp = fits.open(fitsfile)[1].data
        outspec['wave'] = tmp['wave']
        outspec['flux'] = tmp['flux']
        outspec['sig'] = sqrt(1.0/tmp['ivar'])

    elif (instrument == 'MOSFIRE'):
        tmp = fits.open(fitsfile)[1].data
        outspec['wave'] = tmp['wave']
        outspec['flux'] = tmp['flux']
        outspec['sig'] = sqrt(1.0/tmp['ivar'])

    return(outspec)
        
#############################################################################
#############################################################################
#############################################################################
#############################################################################

def loadQsoSpec(obj_id, spectrographs=['XShooter', 'FIRE', 'MOSFIRE', 'HIRES']):

    if (obj_id < 1 or obj_id > 6):
        print("ERROR: Quasar ID must be between 1 and 6")
        return(None)
    
    db = eigerdb.Eigerdb()
    db.getcursor()

    reply = db.query(f"select name from Quasars where id={obj_id}")
    print(f"Object: {reply[0][0]}, quasarid={obj_id}")
    
    ##### Find the spectra that exist in the observations database

    spectra = {}
    
    for spectrograph in spectrographs:

        querystring = """
        SELECT awsbucket,awspath 
        FROM Observations 
        WHERE quasarid={} 
        AND instrument=\'{}\'""".format(obj_id,spectrograph)

        obs = db.query(querystring)

        if (len(obs) > 0):
            print(f"{spectrograph}:")
            local_file = getSpec(obs[0])
            spectrum = parseSpec(local_file,spectrograph)
            spectra[spectrograph] = spectrum
            if (spectrograph == 'XShooter' and False):
                # Grab both the VIS and NIR arms
                spectrum = getSpec(obs[1])
            else:
                continue
            
    db.close()
    
    return(spectra)

