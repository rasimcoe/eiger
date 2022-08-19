import psycopg2
import sys
import boto3
import botocore
import os
from   eiger.Database.SQL import eigerdb
import numpy as np
from astropy.table import Table

def _createAbsComponentTable():

    cmd = """
    CREATE TABLE AbsComponents (
    id SERIAL    PRIMARY KEY,
    quasarid     INTEGER,
    z           FLOAT,
    z_16pct     FLOAT,
    z_84pct     FLOAT,
    b           FLOAT,
    b_16pct     FLOAT,
    b_84pct     FLOAT
    )
    """
    
    edb = eigerdb.Eigerdb()
    edb.getcursor()
    edb.command(cmd,getreply=False)
    edb.close()

def _createAbsIonTable():

    cmd = """
    CREATE TABLE AbsIons (
    id SERIAL     PRIMARY KEY,
    componentid   INT,
    ion           VARCHAR(10),
    N_med         FLOAT,
    N_16pct       FLOAT,
    N_84pct       FLOAT
    )
    """

    edb = eigerdb.Eigerdb()
    edb.getcursor()
    edb.command(cmd,getreply=False)
    edb.close()

def queryAbsLines(quasarid, zmin=0, zmax=10):

    cmd = f"""
    select ion,z,z_16pct,z_84pct,b,b_16pct, b_84pct, n_med, n_16pct, n_84pct 
    from absions JOIN abscomponents on abscomponents.id=componentid 
    where quasarid={quasarid} AND z > {zmin} AND z < {zmax}
    order by z"""

    edb = eigerdb.Eigerdb()
    edb.getcursor()
    reply = edb.query(cmd)
    edb.close()

    t = Table(names=('ion','z','z_16pct','z_84pct','b','b_16pct', 'b_84pct', 'N_med', 'N_16pct', 'N_84pct'),
              dtype=('S2','f4','f4','f4','f4','f4','f4','f4','f4','f4'))

    for row in reply:
        t.add_row(row)
    
    return(t)
