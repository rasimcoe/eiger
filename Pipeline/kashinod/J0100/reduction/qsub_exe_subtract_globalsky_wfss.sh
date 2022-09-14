#!/bin/sh


# Specify node
# node=messier08

# Input directory
input_dir='calibrated_img2_wfss'

# directories for output and log
output_dir=$input_dir'_globalskySubtV3'


output_dir=$input_dir'_globalskySubtV4'  
log_dir=$output_dir'_logs'

# module
module='b'

mkdir $output_dir
mkdir $log_dir

#globalsky_img=globalsky_nrc"$module"long_wfss.fits
#globalsky_img=globalsky_nrc"$module"long_wfss_cal_globalMedSubt_contSubt_kx101_9_Skx21_5.fits  ## V2
#globalsky_img=globalsky_nrc"$module"long_wfss_cal_globalMedSubt_contSubt_kx51_9_Skx21_5.fits   ## V3
globalsky_img=globalsky_nrc"$module"long_wfss_cal_emlineMasked_globalMedSubt_contSubt_kx51_9_Skx21_5.fits  ## V4

echo $globalsky_img' used' >> $log_dir/aaaREADME.txt
cal_fits_files=`ls ./"$input_dir"/jw0124300100[1,2,3,4]_0[2,4]101_000??_nrc"$module"long_cal.fits`
#cal_fits_files=`ls ./"$input_dir"/jw0124300100[1,2,3,4]_0[2,4]101_000??_nrc"$module"long_cal_globalMedSubt_emline.fits`

echo globalsky_img $globalsky_img
echo $cal_fits_files

for fil in $cal_fits_files; do
    
    echo '------------------------------'
    echo 'cal.fits file: '$fil
    basename=`basename $fil .fits`
    echo 'basename: '$basename
    echo 'output_dir: '$output_dir
    job_name=gsky_subt_$basename
    echo 'job name: '$job_name

    sed -e 's%_NODE_%'$node'%g' \
	-e 's%_JOB_NAME_%'$job_name'%g' \
	-e 's%_CAL_FITS_FILE_NAME_%'$fil'%g' \
	-e 's%_OUTPUT_DIR_%'$output_dir'%g' \
	-e 's%_LOG_DIR_%'$log_dir'%g' \
	-e 's%_BASENAME_%'$basename'%g' \
	-e 's%_GLOBAL_SKY_%'$globalsky_img'%g' \
	exe_subtract_globalsky.sh > exe_subtract_globalsky_$basename.sh
	    
    # Execution
    bash exe_subtract_globalsky_$basename.sh

done

