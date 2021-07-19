import psycopg2
import sys
import boto3
import os

class Eigerdb:

    def __init__(self):
        self.connection = False
        self.cursor     = False

    
    def getcursor(self):

        # Note: Users must specify the database username and hostname (or IP)
        # as environment variables in a .bash_profile or similar file. 
        USR=os.getenv('USER')
        ENDPOINT=os.getenv('EIGERDB_SERVER')

        REGION="us-east-2"
        DBNAME="postgres"
        PORT="5432"
        
        # Get the appropriate local credentials
        session = boto3.Session(profile_name='default')
        
        # Start your local client
        client = session.client('rds')
        
        # Generate the temporary password from the AWS server
        token = client.generate_db_auth_token(DBHostname=ENDPOINT, Port=PORT, \
                                              DBUsername=USR, Region=REGION)
        
        try:
            self.connection = psycopg2.connect(host=ENDPOINT, port=PORT, database=DBNAME, \
                                    user=USR, password=token)
            self.cursor = self.connection.cursor()
            print("Connected to the eiger database")
            
        except Exception as e:
            print("Database connection failed: {}".format(e))       


            
    def query(self, querystring, getreply=True):
        self.cursor.execute(querystring)
        if (getreply == True):
            result = self.cursor.fetchall()
        else:
            result = None
        return(result)

    
    def command(self, querystring, getreply=True):
        self.cursor.execute(querystring)
        if (getreply == True):
            result = self.cursor.fetchall()
        else:
            result = None
        self.connection.commit()
        return(result)
    

    def columnNames(self, tablename):
        querystring = \
            "SELECT * from INFORMATION_SCHEMA.COLUMNS where TABLE_NAME=\'"+tablename+ "\'"
        self.cursor.execute(querystring)
        result = self.cursor.fetchall()
        names = [r[3] for r in result]
        return(names)


    def tableNames(self):
        querystring = \
            "SELECT TABLE_NAME from INFORMATION_SCHEMA.TABLES where TABLE_SCHEMA=\'public\' order by table_name"
        self.cursor.execute(querystring)
        result = self.cursor.fetchall()
        names = [r[0] for r in result]
        return(names)


    def close(self):
        self.connection.close()
        
            
