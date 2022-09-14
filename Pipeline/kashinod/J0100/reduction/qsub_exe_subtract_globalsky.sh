#!/bin/sh


# Specify node
# node=messier08

# # Input directory
# input_dir='calibrated_img2_wfss_medSubt'

# # directories for output and log
# output_dir='calibrated_img2_wfss_medSubt2'
# log_dir='calibrated_img2_wfss_medSubt2_logs'

# Input directory
#input_dir='calibrated_img2_wfss_medSubt3'

# directories for output and log
#output_dir='calibrated_img2_wfss_medSubt4'
#log_dir='calibrated_img2_wfss_medSubt4_logs'


# Input directory
input_dir='calibrated_img2_wfss'

# directories for output and log
output_dir='calibrated_img2_wfss_globalskySubtracted'
log_dir=output_dir'_logs'

mkdir $output_dir
mkdir $log_dir

# Set visit and module
visits='1 2 3 4'
modules='a b'
#detectors='1 2 3 4'
detectors='5'

# [IMPORTANT] filter and mask file
#filter='F115W'
#filter='F200W'
#filter='F356W_imaging'
filter='F356W_wfss'

echo 'Filter: '$filter

for visit in $visits; do

    if [ $visit = '1' ]; then 
        if [ "$filter" = 'F115W' ]; then
            pointings=$(seq 1 12)
        elif [ "$filter" = 'F200W' ]; then
            pointings=$(seq 14 25)
        elif [ "$filter" = 'F356W_wfss' ]; then
            pointings=$(seq 1 12)" "$(seq 14 25)
        elif [ "$filter" = 'F356W_imaging' ]; then
            pointings=$(seq 26 28)
        else
            exit
        fi
    elif [ $visit = '2' ]; then 
        if [ "$filter" = 'F115W' ]; then
            pointings=$(seq 29 40)
        elif [ "$filter" = 'F200W' ]; then
            pointings=$(seq 42 53)
        elif [ "$filter" = 'F356W_wfss' ]; then
            pointings=$(seq 29 40)" "$(seq 42 53)
        elif [ "$filter" = 'F356W_imaging' ]; then
	    pointings=$(seq 54 56)
        else
            exit
        fi
    elif [ $visit = '3' ]; then 
	if [ "$filter" = 'F115W' ]; then
            pointings=$(seq 57 68)
        elif [ "$filter" = 'F200W' ]; then
            pointings=$(seq 70 81)
        elif [ "$filter" = 'F356W_wfss' ]; then
            pointings=$(seq 57 68)" "$(seq 70 81)
        elif [ "$filter" = 'F356W_imaging' ]; then
	    pointings=$(seq 82 84)
        else
            exit
        fi
    elif [ $visit = '4' ]; then 
	if [ "$filter" = 'F115W' ]; then
            pointings=$(seq 85 96)
        elif [ "$filter" = 'F200W' ]; then
            pointings=$(seq 98 109)
        elif [ "$filter" = 'F356W_wfss' ]; then
            pointings=$(seq 85 96)" "$(seq 98 109)
        elif [ "$filter" = 'F356W_imaging' ]; then
	    pointings=$(seq 110 112)
        else
            exit
	fi
    else
	echo 'visit must be either 1 2 3 or 4.'
	exit
    fi

    for pointing in $pointings; do
	pointing=`printf "%05d" "$pointing"`
	for module in $modules; do
	    for detector in $detectors; do

		# Global sky image
		if [ "$filter" = 'F115W' ]; then
                    globalsky_img=globalsky_nrc"$module""$detector"_f115w_masked.fits
		elif [ "$filter" = 'F200W' ]; then
                    globalsky_img=globalsky_nrc"$module""$detector"_f200w_masked.fits
		elif [ "$filter" = 'F356W_imaging' ]; then
                    globalsky_img=globalsky_nrc"$module""$detector"_f356w_imaging_masked.fits
		elif [ "$filter" = 'F356W_wfss' ]; then
                    #globalsky_img=globalsky_nrc"$module""$detector"_f356w_wfss.fits
		    globalsky_img=globalsky_medFiltered_withMask_nrc"$module"5.fits
		else
                    exit
		fi


		echo '-------------------------------------------'
		echo 'node: '$node
		echo 'visit pointing module detector: '$visit $pointing $module $detector
	       
		fil=`ls './'$input_dir'/jw0124300100'$visit'_01101_'$pointing'_nrc'$module$detector'_cal.fits'`
		basename=`basename $fil .fits`
		job_name=gsky_subt_$basename
		echo 'cal.fits file:   '$fil
		echo 'globalsky image: '$globalsky_img
		echo 'basename:        '$basename
		echo 'output_dir:      '$output_dir
		echo 'job name:        '$job_name
		
		sed -e 's%_NODE_%'$node'%g' \
		    -e 's%_JOB_NAME_%'$job_name'%g' \
		    -e 's%_CAL_FITS_FILE_NAME_%'$fil'%g' \
		    -e 's%_OUTPUT_DIR_%'$output_dir'%g' \
		    -e 's%_LOG_DIR_%'$log_dir'%g' \
		    -e 's%_BASENAME_%'$basename'%g' \
		    -e 's%_GLOBAL_SKY_%'$globalsky_img'%g' \
		    exe_subtract_globalsky.sh > exe_subtract_globalsky_$basename.sh
	    
	        # Execution
		qsub exe_subtract_globalsky_$basename.sh

	    done
	done
    done
done
