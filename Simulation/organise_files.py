import numpy

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


data_dir='output/'
destination_dir='organised_output/simulation_b100m26z57/'
ensure_dir_exists(destination_dir)


FILTERNAMES=['F115W','F200W','F356W']


#CREATE FOLDERS
for filters in FILTERNAMES:
	if os.path.exists(destination_dir+'IMAGING_%s/'%filters)==False:
		os.mkdir(destination_dir+'IMAGING_%s/'%filters)

if os.path.exists(destination_dir+'GRISM_F356W/')==False:
	os.mkdir(destination_dir+'GRISM_F356W/')

if os.path.exists(destination_dir+'CATALOGS/')==False:
	os.mkdir(destination_dir+'CATALOGS/')

#MOVE CATALOGS
os.system('cp gal_test_grism.cat %s/CATALOGS/'%(destination_dir))
os.system('cp gal_test_im.cat %s/CATALOGS/'%(destination_dir))
os.system('cp ptsrc_test_grism.cat %s/CATALOGS/'%(destination_dir))
os.system('cp ptsrc_test_im.cat %s/CATALOGS/'%(destination_dir))
os.system('cp %s/yamls/test_sed_file.hdf5 %s/CATALOGS/'%(data_dir,destination_dir))


#SELECTING FILELISTS
yaml_files = glob(data_dir+'yamls/jw01243001001*.yaml')

yaml_WFSS_files = []
files_f115w = []
files_f200w = []
files_f356w = []


for f in yaml_files:
    my_dict = yaml.safe_load(open(f))
    if my_dict["Inst"]["mode"]=="wfss":
        yaml_WFSS_files.append(f)
    if my_dict["Inst"]["mode"]=="imaging" and my_dict["Readout"]["filter"]=="F356W":
        files_f356w.append(f)
    if my_dict["Inst"]["mode"]=="imaging" and my_dict["Readout"]["filter"]=="F115W":
        files_f115w.append(f)
    if my_dict["Inst"]["mode"]=="imaging" and my_dict["Readout"]["filter"]=="F200W":
        files_f200w.append(f)


print(len(yaml_WFSS_files))


#NOW MOVING FILES 
FILTER_FILES=[files_f115w,files_f200w,files_f356w]

for q in range(len(FILTERNAMES)):
	THISFILELIST=FILTER_FILES[q]
	for thisfile in THISFILELIST:
		image_name=thisfile[:-5]+'_uncal.fits'
		command='mv %s %sIMAGING_%s/'%(image_name,destination_dir,FILTERNAMES[q])
		command=command.replace('/yamls/','/')
		os.system(command)
		print(command)


for thisfile in yaml_WFSS_files:
		image_name=thisfile[:-5]+'_uncal.fits'
		command='mv %s %sGRISM_F356W/'%(image_name,destination_dir)
		command=command.replace('/yamls/','/')
		print(command)

		os.system(command)
		