from datetime import datetime
import os, sys
from mirage.yaml import yaml_generator

print('MIRAGE_HOME=', os.environ['MIRAGE_HOME'])
print('MIRAGE_DATA=', os.environ['MIRAGE_DATA'])

mirage_dir = os.environ['MIRAGE_HOME']+'/'

# Observation
apt_dir=mirage_dir+'APT_for_mirage/sim/'
pointing_file=apt_dir+'gto1243_2021Oct_J0100.pointing'
xml_file=apt_dir+'gto1243_2021Oct_J0100.xml'

print('DK - ', datetime.now(), ' - input pointing_file: ', pointing_file, flush=True)
print('DK - ', datetime.now(), ' - input xml_file: ', xml_file, flush=True)

# Catalogs
cat_dir = mirage_dir+'catalogs/sourceCatalogs'
cat_name_pts = 'pts_QSO_J0100.txt'
cat_name_gal = 'gal_J0100_z0_15_m999_pa00_r1_wQSO.txt'

print('DK - ', datetime.now(), ' - input quasar catalog: ', os.path.join(cat_dir,cat_name_pts), flush=True)
print('DK - ', datetime.now(), ' - input galaxy catalog: ', os.path.join(cat_dir,cat_name_gal), flush=True)

catalogs = {'J0100+2802': {'nircam':{'point_source':os.path.join(cat_dir,cat_name_pts),
                                     'galaxy': os.path.join(cat_dir,cat_name_gal),
                                     }
                           }
            }

# Background specification
# ADU/sec/pixel
# defalut 'low'
background=None

# Roll angle
roll_angle = None

# Observation dates
# If not set or None, default 2021-10-04
dates = None

# Cosmic ray rates
# If not set or None, default SUNMAX Library with a scaling factor of 1.0
cr = None

# Signal limit for the segmentation map
# This is for WFSS only
# default 0.031 ADU/sec
minimum_signal = 0.031
minimum_signal_units = 'ADU/sec'

# Run the Yaml Generator
print('DK - ', datetime.now(), ' - Start generating yaml files.', flush=True)

yam = yaml_generator.SimInput(input_xml=xml_file, 
                              pointing_file=pointing_file, 
                              catalogs=catalogs, 
                              verbose=True,
                              output_dir='./yaml_files',
                              simdata_output_dir='./simulated_data',
                              cosmic_rays=cr, 
                              background=background, 
                              roll_angle=roll_angle,
                              dates=dates, datatype='linear, raw',
                              segmap_flux_limit=minimum_signal,
                              segmap_flux_limit_units=minimum_signal_units,
                              dateobs_for_background=False,
                              reffile_defaults='crds')

yam.use_linearized_darks = True
yam.create_inputs()

print('DK - ', datetime.now(), ' - Ended generating yaml files.', flush=True)


