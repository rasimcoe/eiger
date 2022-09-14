
import os, sys
import time
from datetime import datetime
from glob import glob


# import pipeline
from jwst.pipeline import Image3Pipeline

# Input asn files
asn_file = sys.argv[1]
print('DK - ', datetime.now(), ': Input asn_file: ', asn_file, flush=True)
if not os.path.isfile(asn_file):
    print('DK - ', datetime.now(), ': not found: ', asn_file, flush=True)
    sys.exit()


output_dir=sys.argv[2]

if not os.path.isdir(output_dir):
    print('mkdir ', output_dir, flush=True)
    os.mkdir(output_dir)
    
# Create instance
pipe = Image3Pipeline()
pipe.output_dir=output_dir
pipe.save_results=True
#pipe.save_calibrated_ramp=True
pipe.skymatch.skip=True
pipe.tweakreg.skip=True

#pipe.resample.save_results=True
#pipe.resample.output_file='reample.fits'
#pipe.resample.output_dir=output_dir
#pipe.outlier_detection.skip=True

filename = asn_file
print('DK - ', datetime.now(), ': Start processing '+filename, flush=True)
pipe.run(filename)
print('DK - ', datetime.now(), ': Ended processing '+filename, flush=True)

#try:
#    pipe.run(filename)
#    print('DK - ', datetime.now(), ': Ended processing '+filename, flush=True)
#except:
#    print('DK - ', datetime.now(), ': Failed processing '+filename, flush=True)
#    raise ValueError


print('DK - ', datetime.now(),': Congratulations, completed!!!', flush=True)

