#!/bin/sh


# Specify node
node=messier11

# Set visit and module
visits='2 3 4'
modules='a b'

# Radius factor
rad_fact='2'

# log directory
output_dir='calibrated_img2_wfss_mask_r2'
log_dir=$output_dir'_logs'

mkdir $output_dir
mkdir $log_dir


for visit in $visits; do

    if [ $visit = '1' ]; then 
	pointings="$(seq 1 12) $(seq 14 25)"
    elif [ $visit = '2' ]; then 
	pointings="$(seq 29 40) $(seq 42 53)"
    elif [ $visit = '3' ]; then 
	pointings="$(seq 57 68) $(seq 70 81)"
    elif [ $visit = '4' ]; then 
	pointings="$(seq 85 96) $(seq 98 109)"
    else
	echo 'visit must be either 1 2 3 or 4.'
	exit
    fi

    for pointing in $pointings; do
	pointing=`printf "%05d" "$pointing"`
	for module in $modules; do

	    mask_fil=sextractor/results/sex_jw01243_nrc_img3_wfss_visit"$visit$module"_i2d_selected.reg
	    mask_fil=`ls $mask_fil`
	    echo '-------------------------------'
	    echo 'node: '$node
	    echo 'visit pointing module detector: '$visit $pointing $module $detector
	    
	    fil=`ls './calibrated_img2_wfss/jw0124300100'$visit'_01101_'$pointing'_nrc'$module'5_cal.fits'`
	    echo 'cal.fits:   '$fil
	    echo 'mask fil:   '$mask_fil
	    basename=`basename $fil .fits`
	    echo 'basename:   '$basename
	    
	    job_name=mask_wfss_$basename
	    echo 'job name:   '$job_name
	    echo 'output dir: '$output_dir
	    echo 'rad_fact:   '$rad_fact

		sed -e 's%_NODE_%'$node'%g' \
		    -e 's%_JOB_NAME_%'$job_name'%g' \
		    -e 's%_CAL_FITS_FILE_NAME_%'$fil'%g' \
		    -e 's%_OUTPUT_DIR_%'$output_dir'%g' \
		    -e 's%_LOG_DIR_%'$log_dir'%g' \
		    -e 's%_BASENAME_%'$basename'%g' \
		    -e 's%_MASK_FIL_%'$mask_fil'%g' \
		    -e 's%_RAD_FACT_%'$rad_fact'%g' \
		    exe_mask_objects_img2cal.sh > exe_mask_objects_img2cal_$basename.sh
	    
	        # Execution
		qsub exe_mask_objects_img2cal_$basename.sh

	done
    done
done
