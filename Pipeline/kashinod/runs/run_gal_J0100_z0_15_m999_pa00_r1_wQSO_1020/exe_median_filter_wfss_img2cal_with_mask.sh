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

kx=_KERNEL_X_
ky=_KERNEL_Y_
kx_gap=_KERNEL_XGAP_

# Execution
echo 'input file: '$fil
echo 'output dir: '$outdir 
echo 'kernel:     '$kx $ky $kx_gap
echo 'mask file:  '$mask_fil
python ../scripts/median_filter_wfss_img2cal_v20220328.py $fil $output_dir $kx $ky $kx_gap --mask $mask_fil


mv _JOB_NAME_.o* $log_dir
mv _JOB_NAME_.e* $log_dir
mv exe_median_filter_wfss_img2cal_with_mask_$basename.sh $log_dir



