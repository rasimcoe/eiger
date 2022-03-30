
import os, sys
import time
from datetime import datetime
from glob import glob
from multiprocessing import Pool

# import pipeline
from jwst.pipeline import Detector1Pipeline

# Input uncal fits files
if len(sys.argv)==3:
    uncal_fits_file = sys.argv[1]
    output_dir = sys.argv[2]
else:
    print('DK - ', datetime.now(), ': [Usage] pipeline_Detector1.py simulated_data/jw...uncal.fits calibrated', flush=True)
    raise ValueError

print('DK - ', datetime.now(), ' : Input uncal.fits file: ', uncal_fits_file, flush=True)

if not os.path.isfile(uncal_fits_file):
    print('DK - ', datetime.now(), ': Not found '+uncal_fits_file, flush=True)
    raise ValueError

print('DK - ', datetime.now(),': Input uncal.fits file found.', flush=True)

# Check if the output directory exits.
if not os.path.isdir(output_dir):
    print('DK - ', datetime.now(),': no directory for output: ', output_dir, flush=True)
    raise ValueError


# Create instance
print('DK - ', datetime.now(),': Creating Detector1Pipeline instance.', flush=True)
pipe = Detector1Pipeline()

pipe.output_dir=output_dir
pipe.save_results=True
pipe.save_calibrated_ramp=True

print('DK - ', datetime.now(), ': Start processing '+uncal_fits_file, flush=True)
try:
    result = pipe.run(uncal_fits_file)
    print('DK - ', datetime.now(), ': Ended processing '+uncal_fits_file, flush=True)
except:
    print('DK - ', datetime.now(), ': Failed processing '+uncal_fits_file, flush=True)
    raise ValueError

print('DK - ', datetime.now(),': Congratulations, completed!!!', flush=True)

