import eigerdb
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from numpy import log10

class MockGalaxy:

    def __init__(self):
        self.id = -1
        self.objname = ' '
        self.qsoid = -1
        self.ra_string = ' '
        self.dec_string = ' '
        self.ra_deg = 0
        self.dec_deg = 0
        self.z = -1
        self.impact_par = -1
        self.f606w = 0
        self.f775w = 0
        self.f814w = 0
        self.f850lp = 0
        self.f115w = 0
        self.f200w = 0
        self.f356w = 0
        self.reff = 0
        self.axis_ratio = 0
        self.pasky = 0
        self.paqso = 0
        self.sersic = 0
        self.m1400 = 0
        self.mstar = 0
        self.sfr = 0
        self.mh = 0
        self.ew_5007 = 0
        self.ew_6563 = 0
        self.ew_4861 = 0

def select(querystring):

    # Example query strings

    # Selects a bunch of fields near the quasar, brighter than 28.75 mag in f775w
    q = """select id, ra_string,dec_string,z,f200w, f356w,impact_par, mstar 
    from MockGalaxies 
    where (impact_par < 10 and f775w < 28.75)"""

    # This is a rectangle box search around the QSO
    q = """select id, ra_deg, dec_deg, z, f356w, impact_par 
    from MockGalaxies 
    where (abs(ra_deg-15.05423)*3600*cos(dec_deg*3.14159/180)<10 
    and abs(dec_deg-28.0405)*3600 < 10)"""

    # This is a circle/cone search around the QSO
    # (This could also be done using the impact_par keyword,
    # but this is a demonstration.
    q = """select id, ra_deg, dec_deg, z, f356w, impact_par 
    from MockGalaxies 
    where (sqrt(((ra_deg-15.05423)*3600*cos(dec_deg*3.14159/180))^2 + 
    ((dec_deg-28.0405)*3600)^2) < 20 and f775w < 27)
    order by impact_par"""

    # Make the search and return results.
    edb = eigerdb.Eigerdb()
    edb.getcursor()
    reply = edb.command(q,getreply=True)
    edb.close()
    return(reply)
        
# Should only be called once, to create the initial table.
def _createMockGalaxiesTable(self):
    create_mocks_table = """
    CREATE TABLE MockGalaxies (
    id SERIAL    PRIMARY KEY,
    objname      VARCHAR ( 15 ),
    qsoid        int,
    ra_string    VARCHAR ( 15 ),
    dec_string   VARCHAR ( 15 ),
    ra_deg       float,
    dec_deg      float,
    z            float,
    impact_par   float,
    f606w        float,
    f775w        float,
    f814w        float,
    f850LP       float,
    f115w        float,
    f200w        float,
    f356w        float,
    reff         float,
    axis_ratio   float,
    pasky        float,
    paqso        float,
    sersic       float,
    m_1400       float,
    mstar        float,
    sfr          float,
    mh           float,
    ew_5007      float,
    ew_6563      float,
    ew_4861      float,
    spec_awsbucket VARCHAR ( 15 ),
    spec_awspath   VARCHAR ( 15 ),
    spec_awsurl    VARCHAR ( 15 )
    )
    """
    
    edb = eigerdb.Eigerdb()
    edb.getcursor()
    edb.command(create_mocks_table,getreply=False)
    edb.close()



# Used to create the catalong database.  Note that the catalog has 87807 entries
# but there was apparently an error in row 87275 that stopped the calculation, so
# there are a few objects that did not make it into the database.

def createMockCatalogFromSimulation(table='/Users/simcoe/Science/eiger/Mirage/input_galaxycatalog_EIGERsim.fits'): 

    hdu = fits.open(table)
    bigtable = hdu[1].data

    # These are the quasar coordinates supplied by Jorryt
    qso_ra  = 15.05423
    qso_dec = 28.0405
    qso_coord = SkyCoord(qso_ra*u.deg, qso_dec*u.deg)

    edb = eigerdb.Eigerdb()
    edb.getcursor()
    
    n_entries = len(bigtable)
    
    for entry in range(n_entries):
        
        objname         = bigtable['ID'][entry]
        qsoid           = 7
        ra_string       = ' '
        dec_string      = ' '
        ra_deg          = bigtable['ra'][entry]
        dec_deg         = bigtable['dec'][entry]
        z               = bigtable['redshift'][entry]
        
        gal_coord       = SkyCoord(ra_deg*u.deg, dec_deg*u.deg)
        tt              = gal_coord.ra.hms
        ra_string       = "{:02}:{:02}:{:5.3f}".format(int(tt.h),int(tt.m),tt.s)
        tt              = gal_coord.dec.dms
        dec_string      = "+{:02}:{:02}:{:5.3f}".format(int(tt.d),int(tt.m),tt.s)
        
        impact_par      = gal_coord.separation(qso_coord).arcsecond
        paqso           = gal_coord.position_angle(qso_coord).to(u.deg).to_value()

        # All fluxes in the table are in nJy
        
        f606w           = nJy_to_AB(bigtable['HST_F606W_fnu'][entry])
        f775w           = nJy_to_AB(bigtable['HST_F775W_fnu'][entry])
        f814w           = nJy_to_AB(bigtable['HST_F814W_fnu'][entry])
        f850lp          = nJy_to_AB(bigtable['HST_F850LP_fnu'][entry])
        f115w           = nJy_to_AB(bigtable['NRC_F115W_fnu'][entry])
        f200w           = nJy_to_AB(bigtable['NRC_F200W_fnu'][entry])
        f356w           = nJy_to_AB(bigtable['NRC_F356W_fnu'][entry])
        reff            = bigtable['Re_maj'][entry]
        axis_ratio      = bigtable['axis_ratio'][entry]
        pasky           = bigtable['position_angle'][entry]
        sersic          = bigtable['sersic_n'][entry]
        m_1400          = bigtable['MUV'][entry]
        mstar           = bigtable['mStar'][entry]
        sfr             = bigtable['sfr_10'][entry]
        ew_5007         = bigtable['O3_5007_EW'][entry]
        ew_6563         = bigtable['HBaA_6563_EW'][entry]
        ew_4861         = bigtable['HBaB_4861_EW'][entry]

        insert_command = """INSERT INTO MockGalaxies
        (objname,qsoid,ra_string, dec_string, ra_deg,dec_deg,z, impact_par,
        f606w,f775w,f814w,f850lp,f115w,f200w,f356w,
        reff,axis_ratio,pasky,paqso,sersic,m_1400,mstar,sfr,
        ew_5007,ew_6563,ew_4861) 
        VALUES ( \'{}\', 7, \'{}\', \'{}\', {:.6f}, {:.6f}, {:.4f}, {:.2f}, 
        {:5.3f}, {:5.3f}, {:5.3f}, {:5.3f}, {:5.3f}, {:5.3f}, {:5.3f}, {:5.3f}, {:5.3f},
        {:5.3f}, {:5.3f}, {:5.3f}, {:5.3f}, {:5.3f}, {:5.3f}, {:5.3f}, {:5.3f}, {:5.3f}
        )""".format(objname,ra_string,dec_string,ra_deg, dec_deg, z, impact_par, \
                    f606w, f775w,f814w,f850lp,f115w,f200w,f356w,reff,axis_ratio, \
                    pasky,paqso,sersic,m_1400,mstar,sfr,ew_5007,ew_6563,ew_4861)
                    
        # print(insert_command)

        edb.command(insert_command, getreply=False)

    edb.close()

def nJy_to_AB(fnu):
    
    return(-2.5 * log10(fnu * 1e-9 * 1e-23) - 48.6)
    
