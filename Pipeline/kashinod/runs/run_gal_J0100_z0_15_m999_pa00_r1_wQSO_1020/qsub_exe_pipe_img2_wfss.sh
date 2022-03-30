#!/bin/sh


# Specify node
node=messier06

# Input file directory
input_dir='calibrated_det1_wfss_hdr_corr'

# Output directory
output_dir='calibrated_img2_wfss'
mkdir $output_dir

log_dir=$output_dir'_logs' 
mkdir $log_dir

# Set visit and module
visits='1 2 3 4'
modules='a b'

for visit in $visits; do

    if [ $visit = '1' ]; then 
	pointings=$(seq 1 12)' '$(seq 14 25)
    elif [ $visit = '2' ]; then 
	pointings=$(seq 29 40)' '$(seq 42 53)
    elif [ $visit = '3' ]; then 
	pointings=$(seq 57 68)' '$(seq 70 81)
    elif [ $visit = '4' ]; then 
	pointings=$(seq 85 96)' '$(seq 98 109)
    else
	echo 'visit must be either 1 2 3 or 4.'
	exit
    fi

    for pointing in $pointings; do
	pointing=`printf "%05d" "$pointing"`
	for module in $modules; do
	    echo '-------------------------------'
	    echo 'node: '$node
	    echo 'visit pointing module: '$visit $pointing $module

	    fil=`ls './'$input_dir'/jw0124300100'$visit'_01101_'$pointing'_nrc'$module'5_rate.fits'`
	    echo 'rate fits file: '$fil

	    basename=`basename $fil .fits`
	    echo 'base name: '$basename
	    echo 'output_dir: '$output_dir
	    job_name=pipe_img2_wfss_$basename
	    echo 'job name: '$job_name

	    sed -e 's%_NODE_%'$node'%g' \
		-e 's%_JOB_NAME_%'$job_name'%g' \
		-e 's%_RATE_FILE_NAME_%'$fil'%g' \
		-e 's%_INPUT_DIR_%'$input_dir'%g' \
		-e 's%_OUTPUT_DIR_%'$output_dir'%g' \
		-e 's%_LOG_DIR_%'$log_dir'%g' \
		-e 's%_BASENAME_%'$basename'%g' \
		exe_pipe_img2.sh > exe_pipe_img2_$basename.sh
	    
	    # Execution
	    qsub exe_pipe_img2_$basename.sh

	done
    done
done
