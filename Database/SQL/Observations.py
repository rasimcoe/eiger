import psycopg2
import sys
import boto3
import os
import eigerdb

class Observation:

    def __init__(self, instrument=None, observatory=None, utdate=None, datafile=None, errfile=None, quasarID=None):
        self.id            = None
        self.quasarid      = quasarID
        self.instrument    = instrument
        self.observatory   = observatory
        self.ut_date       = utdate
        self.exptime       = None
        self.filtername    = None
        self.disperser     = None
        self.awsbucket     = 'gto1243'
        self.awspath       = None
        self.awsurl        = None
        self.awspath_err   = None
        self.awsurl_err    = None
        self.localfile     = datafile
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
        if (True):
            edb = eigerdb.Eigerdb()
            edb.getcursor()
            edb.command(kill_command)
            edb.shutdown()
        else:
            print("ERROR: Please don't do this unless you really know what you are doing")
            print("ERROR: It will destroy the table and all embedded data")
            print("DROP TABLE command not sent")

            
    def addToDatabase(self):

        # Open a connection to the database
        edb = eigerdb.Eigerdb()
        edb.getcursor()
        
        # Check that observatory name is legal
        obscheck = edb.query("select id from observatories where observatory=\'{}\'".format(self.observatory))
        if (len(obscheck) == 0):
            print("ERROR: Observatory is not yet listed in the database, or is formatted inconsistently")
            print("Legal choices are:")
            print(edb.columnNames('observatories'))
            edb.close()
            return False

        # Check that instrument name is legal
        obscheck = edb.query("select id from instruments where name=\'{}\'".format(self.instrument))
        if (len(obscheck) == 0):
            print("ERROR: Instrument is not yet listed in the database, or is formatted inconsistently")
            print("Legal choices are:")
            print(edb.columnNames('instruments'))
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

        ##########  Upload the datafile #############

        try:
            s3_resource.Bucket(self.awsbucket).upload_file(self.localfile,self.awspath)
        except:
            print("Error: file upload aborted")

        ##########  Add this to the SQL table of observations ###########

        query_string = \
        """INSERT INTO observations
        (instrument,observatory,filter,exptime,awsbucket,awspath)
        VALUES (\'{}\',\'{}\',\'{}\',{},\'{}\',\'{}\')
        """.format(self.instrument,self.observatory,self.filter,self.exptime,self.awsbucket,self.awspath)
        edb.command(query_string,getreply=False)


        print("addToDatabase: All done!")
        edb.close()
        
