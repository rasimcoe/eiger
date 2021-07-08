import psycopg2
import sys
import boto3
import os

ENDPOINT="gto1243.clvvehluzo54.us-east-2.rds.amazonaws.com"
PORT="5432"
USR="eigertest"
REGION="us-east-2"
DBNAME="postgres"

# Get the appropriate local credentials
session = boto3.Session(profile_name='default')

# Start your local client
client = session.client('rds')

# Generate the temporary password from the AWS server
token = client.generate_db_auth_token(DBHostname=ENDPOINT, Port=PORT, \
                                      DBUsername=USR, Region=REGION)

try:
    conn = psycopg2.connect(host=ENDPOINT, port=PORT, database=DBNAME, \
                            user=USR, password=token)

    cur   = conn.cursor()
    reply = cur.execute("""SELECT id,name from quasars""")
    query_results = cur.fetchall()
    print(query_results)

except Exception as e:
    print("Database connection failed: {}".format(e))                
                
conn.close()
