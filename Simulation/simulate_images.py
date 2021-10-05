#all these imports are standard for mirage // based on tutorials
from glob import glob
import os
import shutil
import urllib
#These tell mirage the location of calibration data
os.environ["MIRAGE_DATA"] = "/scratch/EIGER/simulate_mirage/MIRAGE_DATA/mirage_data"
os.environ["CRDS_DATA"] = "/scratch/EIGER/simulate_mirage/MIRAGE_DATA/crds_cache"
os.environ["CRDS_PATH"] = "/scratch/EIGER/simulate_mirage/MIRAGE_DATA/crds_cache"
os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"

import pysiaf
import h5py
from astropy.io import ascii as asc
from astropy.io import fits
from matplotlib import cm
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt

from mirage import imaging_simulator
from mirage.catalogs import create_catalog
from mirage.utils.utils import ensure_dir_exists
from mirage.yaml import yaml_generator
import yaml
from mirage.catalogs import catalog_generator
from astropy.io import fits

#these are our own "Eiger" inputs
import eiger_mirage

###VARIABLE INPUT PARAMETERS
pointing_file = '../apt_files/1243_single_backup.pointing'  #Generated with APT and our program ID (1243) - for a single QSO
xml_file = '../apt_files/1243_single_backup.xml' #Generated with APT and our program ID (1243) - for a single QSO
output_dir='output'
ensure_dir_exists(output_dir)


SW=False #Set to True if you want to simulate SW F115W and F200W as well (otherwise it will only simulate F356W images)
CREATE_YAML=False #Set to False if you already created yaml files

target_RA=15.05423  #The RA of the center of the pointing
target_DEC=28.0405  #The DEC of the center of the pointing
target_name='J0100+2802' #Target name
MAGLIM=26. #The magnitude limit in either filter of sources to include in the simulation (this means that sources with magnitudes 31,30.5,29.5 in F115W,F200W and F356W are still included, for example)
FLUX_MULTIPLY=1. #Arteficially brighten sources by this factor - only relevant for grism and ignored for images
redshift_lim_low=5. #Lower limit of the redshift of sources to include 
redshift_lim_high=7. #Lower limit of the redshift of sources to include 


zqso=6.3 #Redshift of the QSO
AB_MAGNITUDES=[17.9,17.0,19.7] #AB magnitudes of the QSO in the filters J,Ks,W1

W1_QSO=19.7 #AB magnitude in the filter closest to JWST F356W - in this case Wise 1 filter -- Only relevant for Grism, not imaging
lamb_W1=33526.00 #in Angstrom  -- Only relevant for Grism, not imaging
average_flamb_QSO=(3.34E4)**-1 * (lamb_W1)**-2 * 10**(-0.4*(W1_QSO-8.9)) #this is in erg/s/cm2/A -- Only relevant for Grism, not imaging

###END OF VARIABLE INPUT PARAMETERS


###CATALOG CREATION - to be used as seed catalog
#Create Pointsource catalog -- in this case simply the QSO
eiger_mirage.create_pointsource_catalog(target_RA,target_DEC,AB_MAGNITUDES,filename='ptsrc_test_im.cat')

#Create galaxy catalog:
eiger_mirage.create_galaxy_catalog('../catalog_data/JADES_SF_mock_r1_v1.2.fits',target_RA,target_DEC,MAGLIM,redshift_lim_low,redshift_lim_high,filename='gal_test_im.cat')

cat_dict = {target_name: {'galaxy':'gal_test_im.cat','point_source':'ptsrc_test_im.cat'}}   
###END OF CATALOG CREATION




###CREATION OF YAML FILES - these ascii text files contain all relevant information that mirage needs to simulate an image. Note that both imaging & grism yaml files are created now
# Create a series of data simulator input yaml files from APT files
yaml_dir = output_dir+'/yamls'
ensure_dir_exists(yaml_dir)


yam = yaml_generator.SimInput(input_xml=xml_file, pointing_file=pointing_file,catalogs=cat_dict,verbose=True, output_dir=yaml_dir, simdata_output_dir=output_dir, segmap_flux_limit=0.01,segmap_flux_limit_units='ADU/sec',datatype='linear,raw') 
#Notes:
#standard segmap_flux_limit is 0.031. A lower limit means that sources are included out to a lower surface brightness. 
#Datatype 'raw' alone is enough for the pipeline. 'linear' won't be used as input for the pipeline but it's useful for visually checking

if CREATE_YAML==True:
	yam.create_inputs()	#If you already created the yaml files you can comment this out
###END CREATION OF YAML FILES




###CREATING SIMULATED IMAGES
yaml_files = glob(yaml_dir+'/jw01243001001*.yaml') #Only going to do visit 1 for now

#Splitting the list by filetype - this is not really necessary but useful for only simulating a certain type of images
F356W_yaml_imaging_files = []
F115W_yaml_imaging_files = []
F200W_yaml_imaging_files = []

for f in yaml_files:
    my_dict = yaml.safe_load(open(f))
    if my_dict["Inst"]["mode"]=="imaging" and my_dict["Readout"]["filter"]=="F356W":
        F356W_yaml_imaging_files.append(f)
    if my_dict["Inst"]["mode"]=="imaging" and my_dict["Readout"]["filter"]=="F115W":
        F115W_yaml_imaging_files.append(f)
    if my_dict["Inst"]["mode"]=="imaging" and my_dict["Readout"]["filter"]=="F200W":
        F200W_yaml_imaging_files.append(f)
                
    
#SIMULATE THE IMAGES, Separately by filter so it's easy to only simulate a subset
for thisimage in F356W_yaml_imaging_files: ###Remove the [:2] - this is just so it only simulates the first two images for simplicity
	img_sim = imaging_simulator.ImgSim()
	img_sim.paramfile = thisimage
	img_sim.create()    

#Removing some output products to save diskapce
os.system('rm %s/*seed*image.fits'%output_dir)
os.system('rm %s/*dark_prep*fits'%output_dir) 

#Same for F115W and F200W

if SW==True:
	for thisimage in F115W_yaml_imaging_files:
		img_sim = imaging_simulator.ImgSim()
		img_sim.paramfile = thisimage
		img_sim.create()
	os.system('rm %s/*seed*image.fits'%output_dir)
	os.system('rm %s/*dark_prep*fits'%output_dir) 

	for thisimage in F200W_yaml_imaging_files:
		img_sim = imaging_simulator.ImgSim()
		img_sim.paramfile = thisimage
		img_sim.create()    
  
	os.system('rm %s/*seed*image.fits'%output_dir)
	os.system('rm %s/*dark_prep*fits'%output_dir) 
