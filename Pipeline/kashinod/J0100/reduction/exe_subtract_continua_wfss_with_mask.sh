#!/bin/sh

job_name=_JOB_NAME_
basename=_BASENAME_
fil=_INP_FITS_FILE_NAME_
maskfil=_MASK_FILE_NAME_
output_dir=_OUTPUT_DIR_
log_dir=_LOG_DIR_

kx=_KERNEL_X_
ky=_KERNEL_Y_
kx_gap=_KERNEL_XGAP_

short_kernel=_SHORT_KERNEL_
kx_s=_SKERNEL_X_
ky_s=_SKERNEL_Y_
kx_gap_s=_SKERNEL_XGAP_

# Execution
echo 'input file: '$fil
echo 'output dir: '$output_dir 

if [ $short_kernel == 'True' ]; then
    nohup python subtract_continua_wfss.py $fil $output_dir $kx $ky $kx_gap --mask $maskfil --short_kernel $kx_s $ky_s $kx_gap_s &> $log_dir/$job_name'.txt' &
else
    nohup python subtract_continua_wfss.py $fil $output_dir $kx $ky $kx_gap --mask $maskfil &> $log_dir/$job_name'.txt' &
fi

