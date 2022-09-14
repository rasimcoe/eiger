#!/bin/bash

# Output directory
output_dir='calibrated_img2photomskip_wfss'
mkdir $output_dir

log_dir=$output_dir'_logs' 
mkdir $log_dir


rate_fits_files=`ls ./calibrated_det1_wfss_hdr_corr/jw0124300100[2,3,4]_0[2,4]101_00???_nrc[a,b]long_rate.fits`  ## Visit1
#rate_fits_files=`ls ./calibrated_det1_wfss_hdr_corr/jw01243001001_0[2,4]101_00???_nrc[a,b]long_rate.fits`  ## Visit1
#rate_fits_files=`ls ./calibrated_det1_wfss_hdr_corr/jw01243001002_0[2,4]101_00???_nrc[a,b]long_rate.fits`  ## Visit2
#rate_fits_files=`ls ./calibrated_det1_wfss_hdr_corr/jw01243001003_0[2,4]101_00???_nrc[a,b]long_rate.fits`  ## Visit3
#rate_fits_files=`ls ./calibrated_det1_wfss_hdr_corr/jw01243001004_0[2,4]101_00???_nrc[a,b]long_rate.fits`  ## Visit4

echo ${rate_fits_files}

for fil in $rate_fits_files; do
    echo '---------------------------'
    echo 'node: '$node
    echo 'input file: '$fil
    echo 'output dir: '$output_dir
    echo 'log dir: '$log_dir

    basename=`basename $fil .fits`
    job_name=pipe_det1_$basename
    echo 'job name: '$job_name

    sed -e 's%_NODE_%'$node'%g' \
        -e 's%_JOB_NAME_%'$job_name'%g' \
        -e 's%_FILE_NAME_%'$fil'%g' \
        -e 's%_OUTPUT_DIR_%'$output_dir'%g' \
        -e 's%_LOG_DIR_%'$log_dir'%g' \
        -e 's%_BASENAME_%'$basename'%g' \
        exe_pipe_img2photomskip.sh > $log_dir/exe_pipe_img2photomskip_$basename.sh

    ## Execution
    bash $log_dir/exe_pipe_img2photomskip_$basename.sh

done





