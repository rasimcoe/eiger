


# directory for output and logs
input_dir='calibrated_img2_wfss'
output_dir='calibrated_img2_wfss_medSubt'

input_dir='calibrated_img2_wfss_globalskySubtracted'
output_dir='calibrated_img2_wfss_globalskySubtracted_medSubt'

log_dir=$output_dir'_logs'
mkdir $output_dir
mkdir $log_dir

# Kernel (the width and height)
# Must be 
kernel_x=51
kernel_y=1
kernel_x_gap=9



#cal_fits_files=`ls ./"$input_dir"/jw01243001001_0[2,4]101_000??_nrc[a,b]long_cal.fits`
#cal_fits_files=`ls ./"$input_dir"/jw01243001002_0[2,4]101_000??_nrc[a,b]long_cal.fits`
#cal_fits_files=`ls ./"$input_dir"/jw01243001004_0[2,4]101_000??_nrc[a,b]long_cal.fits`
cal_fits_files=`ls ./"$input_dir"/jw0124300100[1,2,4]_0[2,4]101_000??_nrc[a,b]long_cal.fits`


echo ${cal_fits_files}

for fil in $cal_fits_files; do

    echo '-------------------------------'
    echo 'visit pointing module detector: '$visit $pointing $module
    
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
    bash exe_median_filter_wfss_img2cal_$basename.sh
    
done

