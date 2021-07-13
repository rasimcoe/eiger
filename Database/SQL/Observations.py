import psycopg2
import sys
import boto3
import os
import eigerdb

class Observation:

    def __init__(self, instrument=None, observatory=None, utdate=None, datafile=None, errfile=None):
        self.id            = None
        self.instrument    = instrument
        self.observatory   = observatory
        self.ut_date       = utdate
        self.awsbucket     = 'gto1243'
        self.awspath       = None
        self.awsurl        = None
        self.awspath_err   = None
        self.awsurl_err    = None
        self.localfile     = datafile
        self.localerr      = errfile
        
    def addToDatabase(self):

        # NB: This is an initial implementation, there needs to be much more
        # error checking and exception handling built in to this function
        
        # These are the remote paths on S3, the local file should be passed in argv
        datapath = 'QuasarSpectra/{}/{}'.format(self.instrument,self.localfile)
        errpath = 'QuasarSpectra/{}/{}'.format(self.instrument,self.localerr)

        print("Local:  {}".format(self.localfile))
        print("Remote: {}".format(datapath))
        
        # Get the appropriate local credentials
        s3_resource = boto3.resource('s3')

        # Upload the datafile
        try:
            s3_resource.Bucket(self.awsbucket).upload_file(self.localfile,datapath)
        except:
            print("Error: file upload aborted")

