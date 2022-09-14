#!/bin/sh

job_name=_JOB_NAME_
basename=_BASENAME_
fil=_INP_FITS_FILE_NAME_
maskfil=_MASK_FILE_NAME_
output_dir=_OUTPUT_DIR_
log_dir=_LOG_DIR_


# Execution
echo 'input file: '$fil
echo 'mask file:  '$maskfil
echo 'output dir: '$output_dir 

nohup python subtract_globalMed_wfss.py $fil $output_dir --mask $maskfil &> $log_dir/$job_name'.txt' &




