#!/bin/sh

## uncal.fits file name
fil=_FILE_NAME_

# output directory
output_dir=_OUTPUT_DIR_

log_dir=_LOG_DIR_


echo 'input file: '$fil
echo 'output dir: '$output_dir 

nohup python pipeline_Image2.py $fil $output_dir &> $log_dir/log__JOB_NAME_'.txt' &




