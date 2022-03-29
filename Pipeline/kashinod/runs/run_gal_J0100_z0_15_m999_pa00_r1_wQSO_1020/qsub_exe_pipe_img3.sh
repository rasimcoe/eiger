#!/bin/sh


# Specify node
node=messier06


# asn file
#asn_file=jw01243_nrc_img3_f356w_medSubt_v3_asn.json
#asn_file=jw01243_nrc_img3_f115w_visit4_medSubt_v3_asn.json
#asn_file=jw01243_nrc_img3_f115w_visit3_medSubt_v3_asn.json
asn_file=jw01243_nrc_img3_f115w_nosubpixeldither_medSubt_v3_asn.json
asn_file=jw01243_nrc_img3_f200w_nosubpixeldither_medSubt_v3_asn.json



# Output directory
output_dir='calibrated_img3_medSubt_v3'

# Source catalog arguments (not yet implimented)
# Default values (kernel_fwhm=2.0, snr_threshold=3.0, npixels=5)
# Deblending can work only on axion
kernel_fwhm=2.0
snr_threshold=3.0
npixels=5
deblend=0


# -------------------------------------
echo '-------------------------------'
echo 'node: '$node
echo 'asn file: '$asn_file
echo 'output_dir: '$output_dir
mkdir $output_dir

log_dir=$output_dir'_logs' 
echo 'log_dir: '$log_dir
mkdir $log_dir


basename=`basename $asn_file .json`
echo 'asn basename: '$basename
echo 'source_catalog arguments::::'
echo '  kernel_fwhm:      '$kernel_fwhm
echo '  snr_threshold:    '$snr_threshold
echo '  npixels:          '$npixels
echo '  deblend(0=F/1=T): '$deblend

job_name=pipe_img3_$basename
echo 'job name: '$job_name

sed -e 's%_NODE_%'$node'%g' \
    -e 's%_JOB_NAME_%'$job_name'%g' \
    -e 's%_ASN_FILE_%'$asn_file'%g' \
    -e 's%_OUTPUT_DIR_%'$output_dir'%g' \
    -e 's%_LOG_DIR_%'$log_dir'%g' \
    -e 's%_BASENAME_%'$basename'%g' \
    -e 's%_KERNEL_FWHM_%'$kernel_fwhm'%g' \
    -e 's%_SNR_THRESHOLD_%'$snr_threshold'%g' \
    -e 's%_NPIXELS_%'$npixels'%g' \
    -e 's%_DEBLEND_%'$deblend'%g' \
    exe_pipe_img3.sh > exe_pipe_img3_$basename.sh
	    
# Execution
qsub exe_pipe_img3_$basename.sh
