from glob import glob
import os
import shutil
import urllib
os.environ["MIRAGE_DATA"] = "/scratch/EIGER/simulate_mirage/MIRAGE_DATA/mirage_data"
os.environ["CRDS_DATA"] = "/scratch/EIGER/simulate_mirage/MIRAGE_DATA/crds_cache"
os.environ["CRDS_PATH"] = "/scratch/EIGER/simulate_mirage/MIRAGE_DATA/crds_cache"
os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"
# Third Party Imports
import pysiaf
#all these imports are standard for mirage // based on tutorials
import h5py
from astropy.io import ascii as asc
from astropy.io import fits
from matplotlib import cm
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt

# Local Imports (from nircam_simulator package)
from mirage import imaging_simulator
from mirage.catalogs import create_catalog
from mirage.utils.utils import ensure_dir_exists
from mirage.yaml import yaml_generator
import yaml
from mirage.catalogs import catalog_generator
from astropy.io import fits
from mirage import wfss_simulator
from mirage.apt import read_apt_xml
import astropy.units as u
from mirage.utils.constants import FLAMBDA_CGS_UNITS, FLAMBDA_MKS_UNITS, FNU_CGS_UNITS 

#these are our own "Eiger" inputs
import eiger_mirage




###VARIABLE INPUT PARAMETERS
pointing_file = '../apt_files/1243_single_backup.pointing'  #Generated with APT and our program ID (1243) - for a single QSO
xml_file = '../apt_files/1243_single_backup.xml' #Generated with APT and our program ID (1243) - for a single QSO
output_dir='output'
yaml_dir='output/yamls'
ensure_dir_exists(output_dir)
ensure_dir_exists(yaml_dir)


SIMULATE_GRISM=True #Set to True in case you want to simulate the grism images at the end of this code. Note it's slightly faster to add more multiprocessing (using multiprocess_grism.sh) after running the catalog creation with this code

target_RA=15.05423  #The RA of the center of the pointing
target_DEC=28.0405  #The DEC of the center of the pointing
target_name='J0100+2802' #Target name
MAGLIM=26. #The magnitude limit in either filter of sources to include in the simulation. For grism data this is the F356W magnitude.
FLUX_MULTIPLY=100. #Arteficially brighten sources by this factor - only relevant for grism and ignored for images
redshift_lim_low=5. #Lower limit of the redshift of sources to include 
redshift_lim_high=7. #Lower limit of the redshift of sources to include 

qso_spec='../catalog_data/RESTFRAME_QSO_TEMPLATE.fits' #directory with the QSO template spectrum

zqso=6.3 #Redshift of the QSO
AB_MAGNITUDES=[17.9,17.0,19.7] #AB magnitudes of the QSO in the filters J,Ks,W1

W1_QSO=19.7 #AB magnitude in the filter closest to JWST F356W - in this case Wise 1 filter -- Only relevant for Grism, not imaging
lamb_W1=33526.00 #in Angstrom  -- Only relevant for Grism, not imaging
average_flamb_QSO=(3.34E4)**-1 * (lamb_W1)**-2 * 10**(-0.4*(W1_QSO-8.9)) #this is in erg/s/cm2/A -- Only relevant for Grism, not imaging

###END OF VARIABLE INPUT PARAMETERS


#FIRST GOING TO CHECK IF THE YAML FILES EXIST
yaml_files = glob(yaml_dir+'/jw*.yaml')

yaml_WFSS_files = []  #select the Grism files
for f in yaml_files:
    my_dict = yaml.safe_load(open(f))
    if my_dict["Inst"]["mode"]=="wfss":
        yaml_WFSS_files.append(f)



if len(yaml_WFSS_files)==0:
	print('--->>>> Going to stop the code because there are no WFSS yaml files. Please create them with simulate_images.py <<<<------')
	stop_read_a_few_lines_above



#CATALOG CREATION
##SPECTRUM OF THE QSO
lambda_qso,flux_qso=eiger_mirage.get_qso_spectrum(qso_spec,zqso,average_flamb_QSO)

#CREATE POINTSOURCE CATALOG
eiger_mirage.create_pointsource_catalog(target_RA,target_DEC,AB_MAGNITUDES,filename='ptsrc_test_im.cat')

#CREATE GALAXY CATALOG
sel_spectral_cats=eiger_mirage.create_galaxy_catalog_grism('../catalog_data/JADES_SF_mock_r1_v1.2.fits',target_RA,target_DEC,MAGLIM,redshift_lim_low,redshift_lim_high,filename='gal_test_im.cat')

cat_dict = {'J0100+2802': {'point_source':'ptsrc_test_im.cat','galaxy':'gal_test_im.cat'}}
sed_file=yaml_dir+'/test_sed_file.hdf5'

#CREATE HDF5 FILE WITH THE SPECTRA OF QSO AND GALAXIES
eiger_mirage.create_catalog_grism_spectra('../catalog_data/JADES_SF_mock_r1_v1.2_spec_5A_30um_z_',FLUX_MULTIPLY,sel_spectral_cats,lambda_qso,flux_qso,sed_file)


#END OF CATALOG CREATION


if SIMULATE_GRISM==True:
###SIMULATING GRISM IMAGES
	for thisimage in yaml_WFSS_files[:1]:
		m = wfss_simulator.WFSSSim(thisimage, override_dark=None, save_dispersed_seed=True,
                           extrapolate_SED=True, disp_seed_filename=None, source_stamps_file=None,
                           SED_file=sed_file)
		m.create()
	


	 
            
