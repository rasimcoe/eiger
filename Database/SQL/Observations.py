import psycopg2
import sys
import boto3
import botocore
import os
from eiger.Database.SQL import eigerdb
import numpy as np

class Observation:

    def __init__(self, instrument=None, observatory=None, utdate=None, \
                 datafile=None, errfile=None, quasarID=None,exptime=None,disperser=None):
        self.id            = None
        self.quasarid      = quasarID
        self.instrument    = instrument
        self.observatory   = observatory
        self.ut_date       = utdate
        self.exptime       = exptime
        self.disperser     = disperser
        self.awsbucket     = 'gto1243'
        self.awspath       = None
        self.revision      = None
        self.localfile     = datafile
        self.localerr      = None
        self.localerr      = errfile

    def _createObservationsTable(self):
        create_observations_table = """
        CREATE TABLE Observations (
        id SERIAL    PRIMARY KEY,
        instrument   VARCHAR ( 15 ),
        observatory  VARCHAR ( 15 ),
        ut_date      VARCHAR ( 15 ),
        exptime      FLOAT,
        awsbucket    VARCHAR ( 15 ),
        awspath      VARCHAR ( 200 ),
        awsurl       VARCHAR ( 200 ),
        awspath_err  VARCHAR ( 200 ),
        awsurl_err   VARCHAR ( 200 )
        )
        """
        edb = eigerdb.Eigerdb()
        edb.getcursor()
        edb.command(create_observations_table)
        edb.shutdown()

    def _createObservatoriesTable(self):
        create_observatories_table = """
        CREATE TABLE Observatories (
        id SERIAL    PRIMARY KEY,
        observatory  VARCHAR ( 15 )
        );
        INSERT INTO observatories (observatory) VALUES(\'Magellan\');
        INSERT INTO observatories (observatory) VALUES(\'VLT\');
        INSERT INTO observatories (observatory) VALUES(\'Keck\');
        INSERT INTO observatories (observatory) VALUES(\'HST\');
        """
        edb = eigerdb.Eigerdb()
        edb.getcursor()
        edb.command(create_observatories_table,getreply=False)
        edb.close()

    def _createInstrumentsTable(self):
        create_instruments_table = """
        CREATE TABLE Instruments (
        id SERIAL    PRIMARY KEY,
        name  VARCHAR ( 15 )
        );
        INSERT INTO instruments (name) VALUES(\'FIRE\');
        INSERT INTO instruments (name) VALUES(\'XShooter\');
        INSERT INTO instruments (name) VALUES(\'MOSFIRE\');
        INSERT INTO instruments (name) VALUES(\'NIRES\');
        INSERT INTO instruments (name) VALUES(\'HIRES\');
        INSERT INTO instruments (name) VALUES(\'MagE\');
        INSERT INTO instruments (name) VALUES(\'NIRSPEC\');
        INSERT INTO instruments (name) VALUES(\'NIRCAM\');
        INSERT INTO instruments (name) VALUES(\'ACS\');
        INSERT INTO instruments (name) VALUES(\'WFC\');
        INSERT INTO instruments (name) VALUES(\'IMACS\');
        INSERT INTO instruments (name) VALUES(\'LDSS3\');
        INSERT INTO instruments (name) VALUES(\'MUSE\');
        INSERT INTO instruments (name) VALUES(\'LLAMAS\');
        """
        edb = eigerdb.Eigerdb()
        edb.getcursor()
        edb.command(create_instruments_table,getreply=False)
        edb.close()

        
        
    def addColumn(self, colname, coltype):
        addcol_query = """
        ALTER TABLE observations
        ADD COLUMN {} {};
        """.format(colname, coltype)

        edb = eigerdb.Eigerdb()
        edb.getcursor()
        edb.command(addcol_query,getreply=False)
        edb.close()

    def dropColumn(self, colname, confirm=False):
        if (confirm):
            delcol_query = """
            ALTER TABLE observations
            DROP COLUMN IF EXISTS {}
            """.format(colname)
            
            edb = eigerdb.Eigerdb()
            edb.getcursor()
            edb.command(delcol_query,getreply=False)
            edb.columnNames('observations')
            edb.close()
        else:
            print("dropColumn requires explicit confirmation since it deletes all data from a table column")

        
    def _destroyObservationsTable(self):
        kill_command = "DROP TABLE IF EXISTS observations"
        if (False):
            edb = eigerdb.Eigerdb()
            edb.getcursor()
            edb.command(kill_command)
            edb.shutdown()
        else:
            print("ERROR: Please don't do this unless you really know what you are doing")
            print("ERROR: It will destroy the table and all embedded data")
            print("DROP TABLE command not sent")

            
    def addToDatabase(self, test=True):

        # Open a connection to the database
        edb = eigerdb.Eigerdb()
        edb.getcursor()
        
        # Check that observatory name is legal
        obscheck = edb.query("select id from observatories where observatory=\'{}\'".format(self.observatory))
        if (len(obscheck) == 0):
            print("ERROR: Observatory is not yet listed in the database, or is formatted inconsistently")
            print("Legal choices are:")
            print(edb.query('select observatory from observatories'))
            edb.close()
            return False

        # Check that instrument name is legal
        obscheck = edb.query("select id from instruments where name=\'{}\'".format(self.instrument))
        if (len(obscheck) == 0):
            print("ERROR: Instrument is not yet listed in the database, or is formatted inconsistently")
            print("Legal choices are:")
            print(edb.query('select name from instruments'))
            edb.close()
            return False
        
        # These are the remote paths on S3, the local file should be passed in argv
        localfile_stripped = self.localfile.split('/')[-1]
        self.awspath = '{}/{}/{}'.format(self.observatory,self.instrument,localfile_stripped)

        # Only need to do this if there is an explicit error file
        # Not all pipelines produce this.
        if (self.localerr != None):
            localerr_stripped = self.localerr.split('/')[-1]
            self.awspath_err = '{}/{}/{}'.format(self.observatory,self.instrument,localerr_stripped)

        print("Local:  {}".format(self.localfile))
        print("Remote: {}".format(self.awspath))
        
        ########## Get the appropriate local credentials ###########

        s3_resource = boto3.resource('s3')

        ########## Check versioning in the SQL database     ###########
        ########## If a previous reduction exists, archive  ############
        ########## it with a new revision number and keep   ############
        ########## it in the database. Latest version has   ############
        ########## revision=='current' always               ############

        query_string = """
        select revision,awspath from Observations 
        where quasarid={} and 
        instrument=\'{}\' and 
        disperser=\'{}\'""".format(self.quasarid,self.instrument,self.disperser)

        # response to this should be a listing of all reductions of this object in the database
        # using the same instrument and disperser
        resp = np.array(edb.command(query_string,getreply=True))
        if (len(resp) != 0):
            revs = resp[:,0]
            if (len(revs) == 1):
                old_awspath = resp[0,1]
                maxrev      = 0
                newrev      = 1
            else:
                maxrev = max(revs[revs!='current'])
                old_awspath = resp[revs.index('current'),1]                
                newrev = maxrev+1

            new_awspath = f"{old_awspath[:-5].split('_rev')[0]}_rev{newrev}.fits"
                
            print(f"Renaming prior reduction: {old_awspath}-->{new_awspath}")
            cmd = f"update Observations set awspath=\'{new_awspath}\',revision=\'{newrev}\' where awspath=\'{old_awspath}\'"
            if (test==False):
                # This statement changes the SQL database entry awspath to reflect downrev of the N-1 version
                edb.command(cmd,getreply=False)
                # The next 2 statements change the actual filename of the N-1 version on S3
                # (on S3, you can't rename a file, need to copy to a new file and then delete the original)
                s3_resource.Object(self.awsbucket,new_awspath).\
                    copy_from(CopySource={'Bucket':self.awsbucket,'Key':old_awspath})
                s3_resource.Object(self.awsbucket,old_awspath).delete()
            
        ########## Check to see if a file with this name already exists in the database
        ########## If so, then rename the old file and archive with a revision number.
        
        try:
            s3_resource.Object(self.awsbucket,old_awspath).load()
            print(f"File {old_awspath} still exists but should have been renamed, something has gone wrong.")
            return()
        except botocore.exceptions.ClientError as e:
            if (e.response['Error']['Code'] == "404"):
                # print("File does not exist")
                print("OK to proceed")
            else:
                print("Something has gone wrong accessing S3")
                return(False)
                
        ##########  Upload the datafile to AWS #############

        if (test == False):
            try:
                s3_resource.Bucket(self.awsbucket).upload_file(self.localfile,self.awspath)
            except:
                print("Error: file upload aborted")
                
        ########## Build the database entry string according to how many attributes are known ########

        obsfields = {}
        if (self.observatory != None):
            obsfields['observatory'] = self.observatory
        if (self.instrument != None):
            obsfields['instrument'] = self.instrument
        if (self.disperser != None):
            obsfields['disperser'] = self.disperser
        if (self.exptime != None):
            obsfields['exptime'] = self.exptime
        if (self.awsbucket != None):
            obsfields['awsbucket'] = self.awsbucket
        if (self.awspath != None):
            obsfields['awspath'] = self.awspath
        if (self.ut_date != None):
            obsfields['ut_date'] = self.ut_date
        if (self.quasarid != None):
            obsfields['quasarid'] = self.quasarid
        if (self.revision == None):
            obsfields['revision'] = 'current'

        query_string = "INSERT INTO observations ("
        nfields = len(obsfields)
        i = 0
        for keyname in list(obsfields):
            query_string += keyname
            if (i < nfields-1):
                query_string += ','
                i+=1
            else:
                query_string += ') VALUES ('

        i=0
        for keyval in list(obsfields.values()):

            # Strings need quotes in the query, floats and ints do not.
            if (isinstance(keyval,str)):
                query_string += '\''+keyval+'\''
            else:
                query_string += '{}'.format(keyval)
                
            if (i < nfields-1):
                query_string += ','
                i+=1
            else:
                query_string += ')'
                            
        ##########  Add this to the SQL table of observations ###########

        if (test == False):
            edb.command(query_string,getreply=False)
        else:
            print(query_string)
            
        print("addToDatabase: All done!")
        edb.close()
        
