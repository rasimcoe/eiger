#!/bin/sh

#$ -S /bin/sh
#$ -cwd
#$ -V

# # specifing which node
#$ -q all.q@_NODE_

# # number of using cores
#$ -pe openmpi 1

# # name of job
#$ -N _JOB_NAME_

job_name=_JOB_NAME_
basename=_BASENAME_
fil=_CAL_FITS_FILE_NAME_
output_dir=_OUTPUT_DIR_
log_dir=_LOG_DIR_
global_sky=_GLOBAL_SKY_


echo 'input file: '$fil
echo 'global sky: '$global_sky
echo 'output dir: '$output_dir

python ../scripts/subtract_globalsky.py $fil $output_dir $global_sky


mv exe_subtract_globalsky_$basename.sh $log_dir
mv $job_name.e* $log_dir
mv $job_name.o* $log_dir


