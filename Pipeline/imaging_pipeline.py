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

from multiprocessing import Pool
import jwst

from jwst.pipeline import Detector1Pipeline
from jwst.pipeline import Image2Pipeline
from jwst.pipeline import Image3Pipeline

def run_pipe1(filename):
	pipe = Detector1Pipeline()
	pipe.output_dir=output_dir
	pipe.save_results=True
	pipe.save_calibrated_ramp=True
	pipe.run(filename)


def run_pipe2(filename):
	pipe = Image2Pipeline()
	pipe.output_dir=output_dir
	pipe.save_results=True
	pipe.save_bsub=True
	pipe.run(filename)

simulation_name='simulation_b100m26z57'


for FILTER in ['356']:#,'200','356']:

	output_dir='reduced/%s/imaging_F%sW/'%(simulation_name,FILTER)
	ensure_dir_exists(output_dir)


	for visit in [1]:
		fits_files=glob('../simulate_mirage/organised_output/%s/IMAGING_F%sW/jw0124300100%s*uncal.fits'%(simulation_name,FILTER,visit))
		print(fits_files)

		n_procs = 1 # number of cores available

		#Step 1:
		with Pool(n_procs) as pool:
			pool.map(run_pipe1,fits_files)

	   	#Step 2:
		fits_files=glob('%s/*rate.fits'%output_dir)
		with Pool(n_procs) as pool:
			pool.map(run_pipe2,fits_files)



		#Step 3:
		fits_files=glob('%s/jw0124300100%s*cal.fits'%(output_dir,visit))

		print(fits_files,len(fits_files))
		pipe = Image3Pipeline()
		pipe.output_dir=output_dir
		pipe.save_results=True
		pipe.run(fits_files)




