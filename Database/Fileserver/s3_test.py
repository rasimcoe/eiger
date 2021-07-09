import sys
import boto3
import os

# Get the appropriate local credentials
s3_resource = boto3.client('s3')

try:
    reply = s3_resource.list_buckets()
    # reply = s3_resource.list_objects(Bucket='gto1243')
    print(reply['Buckets'])

except Exception as e:
    print("Fileserver connection failed: {}".format(e))                
                

