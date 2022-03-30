#!/bin/sh


# Specify node
node=messier04
# directory for output and logs
output_dir='calibrated_img2_wfss_medSubt'
log_dir='calibrated_img2_wfss_medSubt_logs'
mkdir $output_dir
mkdir $log_dir

# Kernel (the width and height)
# Must be 
kernel_x=51
kernel_y=1
kernel_x_gap=9

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
	    echo 'visit pointing module detector: '$visit $pointing $module
	    
	    fil=`ls './calibrated_img2_wfss/jw0124300100'$visit'_01101_'$pointing'_nrc'$module'5_cal.fits'`
	    echo 'cal.fits file: '$fil
	    echo 'kernel (x, y, x_gap): '$kernel_x $kernel_y $kernel_x_gap

	    basename=`basename $fil .fits`
	    echo 'basename: '$basename

	    job_name=med_subt_$basename
	    echo 'job name: '$job_name
	    
	    sed -e 's%_NODE_%'$node'%g' \
		-e 's%_JOB_NAME_%'$job_name'%g' \
		-e 's%_CAL_FITS_FILE_NAME_%'$fil'%g' \
		-e 's%_OUTPUT_DIR_%'$output_dir'%g' \
		-e 's%_LOG_DIR_%'$log_dir'%g' \
		-e 's%_BASENAME_%'$basename'%g' \
		-e 's%_KERNEL_X_%'$kernel_x'%g' \
		-e 's%_KERNEL_Y_%'$kernel_y'%g' \
		-e 's%_KERNEL_XGAP_%'$kernel_x_gap'%g' \
		exe_median_filter_wfss_img2cal.sh > exe_median_filter_wfss_img2cal_$basename.sh
	    
	    # Execution
	    qsub exe_median_filter_wfss_img2cal_$basename.sh

	done
    done
done
