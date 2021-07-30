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

import astropy.units as u
from mirage.utils.constants import FLAMBDA_CGS_UNITS, FLAMBDA_MKS_UNITS, FNU_CGS_UNITS 
import sys
import numpy

spectra_catalog_name='output/yamls/test_sed_file.hdf5'


visit=sys.argv[1]
module=sys.argv[2]
quartile=int(sys.argv[3])
print('This core is going to do visit, module, quartile',visit,module,quartile)

#Before running this code, make sure that you have run simulate_grism_catalog.py
if os.path.exists(spectra_catalog_name)==False:
	print('--->>>> Going to stop the code because there are is no hdf5 catalog. Please create it with simulate_grism_catalog.py <<<<------')
	stop_read_a_few_lines_above


#We need to make a copy of the hdf5 file for this process specifically as it can only be read by one core simultaneously

if os.path.isfile('output/yamls/test_sed_file_v%s_m%s_q%s.hdf5'%(visit,module,quartile)) == False:
	os.system('cp output/yamls/test_sed_file.hdf5 ./output/yamls/test_sed_file_v%s_m%s_q%s.hdf5'%(visit,module,quartile))



out_dir='output'


yaml_files = glob('output/yamls/jw0124300100%s*nrc%s5.yaml'%(visit,module)) #can be automated more with the proposal ID and target number

splitted=numpy.array_split(yaml_files,4)



yaml_WFSS_files = []

for f in splitted[quartile]:
    my_dict = yaml.safe_load(open(f))
    if my_dict["Inst"]["mode"]=="wfss":
        yaml_WFSS_files.append(f)

print('Im going to do %s files'%len(yaml_WFSS_files))

for thisfile in yaml_WFSS_files:
		print('Running',thisfile)
		m = wfss_simulator.WFSSSim(thisfile, override_dark=None, save_dispersed_seed=True,
                           extrapolate_SED=True, disp_seed_filename=None, source_stamps_file=None,
                           SED_file='output/yamls/test_sed_file_v%s_m%s_q%s.hdf5'%(visit,module,quartile))
		m.create()	
	
#Removing the hdf5 copy to save space	
os.system('rm output/yamls/test_sed_file_v%s_m%s_q%s.hdf5'%(visit,module,quartile))
        
        
