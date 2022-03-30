#!/bin/sh


# Specify node
node=messier07

# Output directory
output_dir='calibrated_det1'
mkdir $output_dir

log_dir=$output_dir'_logs' 
mkdir $log_dir

visits='4'
modules='a b'
#detectors='1 2 3 4'
#detectors='5_image' # for direct imaging
detectors='5'       # for WFSS



for visit in $visits; do

    if [ $visit = '1' ]; then 
	if [ "$detectors" = '5_image' ]; then
	    pointings=$(seq 26 28)
	else
	    pointings=$(seq 1 12)' '$(seq 14 25)
	fi
    elif [ $visit = '2' ]; then 
	if [ "$detectors" = '5_image' ]; then
	    pointings=$(seq 54 56)
	else	    
	    pointings=$(seq 29 40)' '$(seq 42 53)
	fi
    elif [ $visit = '3' ]; then
	if [ "$detectors" = '5_image' ]; then
	    pointings=$(seq 82 84)
	else	    
	    pointings=$(seq 57 68)' '$(seq 70 81)
	fi
    elif [ $visit = '4' ]; then 
	if [ "$detectors" = '5_image' ]; then
	    pointings=$(seq 110 112)
	else	    
	    pointings=$(seq 85 96)' '$(seq 98 109)
	fi
    else
	echo 'visit must be either 1 2 3 or 4.'
	exit
    fi

    for pointing in $pointings; do
	pointing=`printf "%05d" "$pointing"`
	for module in $modules; do
	    for detector in $detectors; do
		echo '---------------------------'
		echo 'node: '$node
		echo 'input file: '$fil
		echo 'output dir: '$output_dir
		echo 'log dir: '$log_dir

		if [ "$detector" = '5_image' ]; then
		    detector=5
		fi
		fil=`ls './simulated_data/jw0124300100'$visit'_01101_'$pointing'_nrc'$module$detector'_uncal.fits'`
		echo 'uncal fits file: '$fil

		basename=`basename $fil .fits`
		echo 'base name: '$basename

		job_name=pipe_det1_$basename
		echo 'job name: '$job_name

		sed -e 's%_NODE_%'$node'%g' \
		    -e 's%_JOB_NAME_%'$job_name'%g' \
		    -e 's%_FILE_NAME_%'$fil'%g' \
		    -e 's%_OUTPUT_DIR_%'$output_dir'%g' \
		    -e 's%_LOG_DIR_%'$log_dir'%g' \
		    -e 's%_BASENAME_%'$basename'%g' \
		    exe_pipe_det1.sh > exe_pipe_det1_$basename.sh

                # Execution
		qsub exe_pipe_det1_$basename.sh

	    done
	done
    done
done



