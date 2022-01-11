from eiger.Database.SQL import eigerdb
import boto3
from astropy.io import fits


def loadSpec(obj_id):

    if (obj_id < 1 or obj_id > 6):
        print("ERROR: Quasar ID must be between 1 and 5")
        return(None)
    
    db = eigerdb.Eigerdb()
    db.getcursor()

    reply = db.query(f"select name from Quasars where id={obj_id}")
    print(f"Object: {reply[0][0]}, quasarid={obj_id}")
    
    ##### Find the spectra that exist

    # XShooter
    xsh = db.query(f"select awsbucket,awspath FROM Observations WHERE quasarid={obj_id} AND instrument=\'XShooter\'")
    if (len(xsh) > 0):
        print("XShooter: FOUND")
    else:
        print("XShooter: NOT FOUND")
    
    # FIRE
    fire = db.query(f"select awsbucket,awspath FROM Observations WHERE quasarid={obj_id} AND instrument=\'FIRE\'") 
    if (len(fire) > 0):
        print("FIRE:     FOUND")
    else:
        print("FIRE:     NOT FOUND")
   
    # MOSFIRE
    mosfire = db.query(f"select awsbucket,awspath FROM Observations WHERE quasarid={obj_id} AND instrument=\'MOSFIRE\'")
    if (len(mosfire) > 0):
        print("MOSFIRE:  FOUND")
    else:
        print("MOSFIRE:  NOT FOUND")


    
    # HIRES
    hires = db.query(f"select awsbucket,awspath FROM Observations WHERE quasarid={obj_id} AND instrument=\'HIRES\'")
    if (len(hires) > 0):
        print("HIRES:    FOUND")
    else:
        print("HIRES:    NOT FOUND")
    
    ######### Download and package the quasar spectrum files #########

    db.close()

    # Get the appropriate local credentials
    s3_client   = boto3.client('s3')
    s3_resource = boto3.resource('s3')
    
    try:
        reply = s3_resource.Bucket(fire[0][0]).download_file(fire[0][1],'tmp.fits')
    except:
        print("ERROR: AWS file could not be retrieved")
        
    firespec = (fits.open("tmp.fits"))[1].data[0]
    
    return(firespec)

