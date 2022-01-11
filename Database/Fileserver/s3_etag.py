from eiger.Database.SQL import eigerdb
import boto3, io
from hashlib import md5

####################################################################
#
# Each file on S3 has a unique MD5 hash associated with its bits.
# This subroutine pulls the etag from the remote server associated
# with the given file (for comparison with the local value).
# For files smaller than 5Mb, the etag is the MD5 hash.  For large files
# it may be slightly different.

def etag_remotehash(bucket, awspath):

    s3_client   = boto3.client('s3')
    # s3_resource = boto3.resource('s3')

    try:
        reply = s3_client.head_object(Bucket=bucket,Key=awspath)
        etag  = reply['ResponseMetadata']['HTTPHeaders']['etag']
        return(etag[1:-1])
    except:
        print("ERROR: File not found on AWS")
        return('')

####################################################################
#
# This routine takes in a file on local disk and generates an MD5 hash,
# so that we can compare with the cloud-based version of the file. 
# NB for large files there may be issues here, it's more complex than a
# simple MD5 hash in that case.
    
def etag_localhash(filename):

    try:
        etag = md5(io.open(filename,'rb').read()).hexdigest()
        return(etag)
    except:
        print(f"ERROR: file {filename} not found or could not hash")
        return('')

