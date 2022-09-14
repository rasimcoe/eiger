
import os, sys
import time
from datetime import datetime
from glob import glob
from multiprocessing import Pool

# import pipeline
from jwst.pipeline import Image2Pipeline

# Input rate fits files
if len(sys.argv)==3:
    rate_fits_file = sys.argv[1]
    output_dir = sys.argv[2]
else:
    #print('DK - ', datetime.now(), ': [Usage] pipeline_Image2.py calibrated/jw...rate.fits calibrated', flush=True)
    raise ValueError

print('DK - ', datetime.now(), ' : Input rate.fits file: ', rate_fits_file, flush=True)

if not os.path.isfile(rate_fits_file):
    print('DK - ', datetime.now(), ': Not found '+rate_fits_file, flush=True)
    raise ValueError

print('DK - ', datetime.now(),': Input rate.fits file found.', flush=True)

# Check if the output directory exits.
if not os.path.isdir(output_dir):
    print('DK - ', datetime.now(),': no directory for output: ', output_dir, flush=True)
    raise ValueError


# Create instance
print('DK - ', datetime.now(),': Creating Image2Pipeline instance.', flush=True)
pipe = Image2Pipeline()

pipe.output_dir=output_dir
pipe.save_results=True
pipe.save_calibrated_ramp=True

print('DK - ', datetime.now(), ': Start processing '+rate_fits_file, flush=True)
try:
    result = pipe.run(rate_fits_file)
    print('DK - ', datetime.now(), ': Ended processing '+rate_fits_file, flush=True)
except:
    print('DK - ', datetime.now(), ': Failed processing '+rate_fits_file, flush=True)
    raise ValueError

print('DK - ', datetime.now(),': Congratulations, completed!!!', flush=True)


