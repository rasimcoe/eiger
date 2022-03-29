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


kx=_KERNEL_X_
ky=_KERNEL_Y_
kx_gap=_KERNEL_XGAP_

# Execution
echo 'input file: '$fil
echo 'output dir: '$outdir 
echo 'kernel: '$kx $ky
#python ../scripts/median_filter_wfss_img2cal_fits_v20211105.py $fil $output_dir $kx $ky $kx_gpa
python ../scripts/median_filter_wfss_img2cal_v20220328.py $fil $output_dir $kx $ky $kx_gap


mv _JOB_NAME_.o* $log_dir
mv _JOB_NAME_.e* $log_dir
mv exe_median_filter_wfss_img2cal_$basename.sh $log_dir



