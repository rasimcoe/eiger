import sys
import boto3
import os

# Get the appropriate local credentials
s3_client   = boto3.client('s3')
s3_resource = boto3.resource('s3')

try:
    reply = s3_client.list_buckets()
    print(reply['Buckets'])

    remote_filename = 'QuasarSpectra/FIRE/ULAS1120_F.fits'
    local_filename  = 'ULAS1120_F.fits'
    
    reply = s3_resource.Bucket('gto1243').\
        download_file(remote_filename, local_filename)
    
except Exception as e:
    print("Fileserver connection failed: {}".format(e))                
                

