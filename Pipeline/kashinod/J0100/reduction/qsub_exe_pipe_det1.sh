#!/bin/bash


# Output directory
output_dir='calibrated_det1_imaging'
mkdir $output_dir

log_dir=$output_dir'_logs' 
mkdir $log_dir


#uncal_fits_files=`ls ../data/uncal_F356W_GRISM/jw01243001001_02101_000??_nrcalong_uncal.fits`  # WFSS Visit 1 ModA with F115W
#uncal_fits_files=`ls ../data/uncal_F356W_GRISM/jw01243001001_02101_000??_nrcblong_uncal.fits`  # WFSS Visit 1 ModB with F115W
#uncal_fits_files=`ls ../data/uncal_F356W_GRISM/jw01243001001_04101_000??_nrc[a,b]long_uncal.fits`  # WFSS Visit 1 ModA,B with F200W
#uncal_fits_files=`ls ../data/uncal_F356W_GRISM/jw01243001002_0[2,4]101_000??_nrc[a,b]long_uncal.fits`  # WFSS Visit 2
#uncal_fits_files=`ls ../data/uncal_F356W_GRISM/jw01243001003_0[2,4]101_000??_nrc[a,b]long_uncal.fits`  # WFSS Visit 3
#uncal_fits_files=`ls ../data/uncal_F356W_GRISM/jw01243001004_0[2,4]101_000??_nrc[a,b]long_uncal.fits`  # WFSS Visit 4

# Imaging
uncal_fits_files=`ls ../data/uncal_F356W/jw0124300100[2,3,4]_0[5,7]101_0000?_nrc[a,b]long_uncal.fits`  # WFSS Visit 3
#uncal_fits_files=`ls ../data/uncal_F356W/jw0124300100[1,2,3,4]_0[5,7]101_0000[1,2]_nrc[a,b]long_uncal.fits`  # WFSS Visit 3





echo ${uncal_fits_files}

for fil in $uncal_fits_files; do
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
        exe_pipe_det1.sh > $log_dir/exe_pipe_det1_$basename.sh

    ## Execution
    bash $log_dir/exe_pipe_det1_$basename.sh

done





