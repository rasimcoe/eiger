

# Kernel (the width and height)
# Must be 
kernel_x=51     ## default 51
kernel_y=1      ## default 1
kernel_x_gap=9  ## default 9

# Short kernel?
short_kernel=True  ## default True
kernel_x_s=21      ## default 21
kernel_y_s=1       ## default 1
kernel_x_gap_s=5   ## default 5


# directory for output and logs
#input_dir='calibrated_det1_wfss_hdr_corr_globalMedSubt'
#input_dir='calibrated_img2_wfss_globalMedSubt'
#input_dir='calibrated_img2_wfss_globalskySubtV3_globalMedSubt'
input_dir='calibrated_img2_wfss_globalskySubtV4_emlineMasked_globalMedSubt'

if [ $short_kernel == 'True' ]; then
    output_dir=$input_dir'_contSubt_kx'$kernel_x'_'$kernel_x_gap'_Skx'$kernel_x_s'_'$kernel_x_gap_s
else
    output_dir=$input_dir'_contSubt_kx'$kernel_x'_'$kernel_x_gap
fi

log_dir=$output_dir'_logs'
mkdir $output_dir
mkdir $log_dir

## Don't delete, this is for all grism frames.
inp_fits_files=`ls ./"$input_dir"/jw0124300100[1,2,3,4]_0[2,4]101_000??_nrc[a,b]long_cal_globalMedSubt.fits`
inp_fits_files=`ls ./"$input_dir"/jw01243001004_0[2,4]101_000??_nrc[a,b]long_cal_globalMedSubt.fits`


echo ${inp_fits_files}

for fil in $inp_fits_files; do

    echo '-------------------------------'
    echo 'Input .fits file: '$fil

    basename=`basename $fil .fits`
    echo 'basename: '$basename
    
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
	-e 's%_KERNEL_X_%'$kernel_x'%g' \
        -e 's%_KERNEL_Y_%'$kernel_y'%g' \
        -e 's%_KERNEL_XGAP_%'$kernel_x_gap'%g' \
	-e 's%_SHORT_KERNEL_%'$short_kernel'%g' \
	-e 's%_SKERNEL_X_%'$kernel_x_s'%g' \
	-e 's%_SKERNEL_Y_%'$kernel_y_s'%g' \
	-e 's%_SKERNEL_XGAP_%'$kernel_x_gap_s'%g' \
	exe_subtract_continua_wfss_with_mask.sh > $log_dir/exe_subtract_continua_wfss_with_mask_$basename.sh
    
    # Execution
    bash $log_dir/exe_subtract_continua_wfss_with_mask_$basename.sh
    
done

