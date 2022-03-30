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
mask_fil=_MASK_FIL_
rad_fact=_RAD_FACT_
echo 'input file: '$fil
echo 'output dir: '$output_dir

python ../scripts/mask_objects_img2cal.py $fil $output_dir $mask_fil $rad_fact

mv exe_mask_objects_img2cal_$basename.sh $log_dir
mv $job_name.e* $log_dir
mv $job_name.o* $log_dir


