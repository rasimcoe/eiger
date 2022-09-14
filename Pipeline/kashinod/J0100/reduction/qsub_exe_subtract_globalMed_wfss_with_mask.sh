


# directory for output and logs
#input_dir='calibrated_img2_wfss'
input_dir='calibrated_img2_wfss_globalskySubtV4'

output_dir=$input_dir'_emlineMasked_globalMedSubt'

log_dir=$output_dir'_logs'
mkdir $output_dir
mkdir $log_dir

### rate.fits
#inp_fits_files=`ls ./"$input_dir"/jw0124300100[1,2,3,4]_0[2,4]101_000??_nrc[a,b]long_rate.fits`
#inp_fits_files=`ls ./"$input_dir"/jw01243001003_0[2,4]101_000??_nrc[a,b]long_rate.fits`

### cal.fits
inp_fits_files=`ls ./"$input_dir"/jw01243001001_02101_00001_nrc[a,b]long_cal.fits`
#### Don't delete, this is for all grism frames
inp_fits_files=`ls ./"$input_dir"/jw0124300100[1,2,3,4]_0[2,4]101_000??_nrc[a,b]long_cal.fits`

echo ${inp_fits_files}

for fil in $inp_fits_files; do
    
    echo '-------------------------------'
    echo 'Input .fits file: '$fil
    
    basename=`basename $fil .fits`
    echo 'basename: '$basename

    #maskfil=${fil/cal.fits/mask.fits}
    maskfil='calibrated_img2_wfss/'${basename:0:35}'mask.fits'
    echo 'mask fil: '$maskfil
    
    job_name=med_subt_$basename
    echo 'job name: '$job_name
    
    sed -e 's%_NODE_%'$node'%g' \
	-e 's%_JOB_NAME_%'$job_name'%g' \
	-e 's%_INP_FITS_FILE_NAME_%'$fil'%g' \
	-e 's%_MASK_FILE_NAME_%'$maskfil'%g' \
	-e 's%_OUTPUT_DIR_%'$output_dir'%g' \
	-e 's%_LOG_DIR_%'$log_dir'%g' \
	-e 's%_BASENAME_%'$basename'%g' \
	exe_subtract_globalMed_wfss_with_mask.sh > $log_dir/exe_subtract_globalMed_wfss_with_mask_$basename.sh
    
    # Execution
    bash $log_dir/exe_subtract_globalMed_wfss_with_mask_$basename.sh
    
done

