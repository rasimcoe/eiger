#!/bin/sh


# Specify node
node=messier11


# asn file
asn_file=jw01243_nrc_img3_wfss_visit1a_medFiltered_withMask_globalskySubtracted_asn.json
asn_file=jw01243_nrc_img3_wfss_visit1b_medFiltered_withMask_globalskySubtracted_asn.json
asn_file=jw01243_nrc_img3_wfss_visit2a_medFiltered_withMask_globalskySubtracted_asn.json
asn_file=jw01243_nrc_img3_wfss_visit2b_medFiltered_withMask_globalskySubtracted_asn.json
asn_file=jw01243_nrc_img3_wfss_visit3a_medFiltered_withMask_globalskySubtracted_asn.json
asn_file=jw01243_nrc_img3_wfss_visit3b_medFiltered_withMask_globalskySubtracted_asn.json
asn_file=jw01243_nrc_img3_wfss_visit4a_medFiltered_withMask_globalskySubtracted_asn.json
asn_file=jw01243_nrc_img3_wfss_visit4b_medFiltered_withMask_globalskySubtracted_asn.json

# Output directory
output_dir='calibrated_img3_wfss_medFiltered_withMask_globalskySubtracted'

# -------------------------------------
basename=`basename $asn_file .json`
job_name=pipe_img3_wfss_$basename
mkdir $output_dir
log_dir=$output_dir'_logs' 
mkdir $log_dir

echo '-------------------------------'
echo 'node:         '$node
echo 'asn file:     '$asn_file
echo 'output_dir:   '$output_dir
echo 'log_dir:      '$log_dir
echo 'asn basename: '$basename
echo 'job name:     '$job_name

sed -e 's%_NODE_%'$node'%g' \
    -e 's%_JOB_NAME_%'$job_name'%g' \
    -e 's%_ASN_FILE_%'$asn_file'%g' \
    -e 's%_OUTPUT_DIR_%'$output_dir'%g' \
    -e 's%_LOG_DIR_%'$log_dir'%g' \
    -e 's%_BASENAME_%'$basename'%g' \
    exe_pipe_img3_wfss.sh > exe_pipe_img3_wfss_$basename.sh
	    
# Execution
qsub exe_pipe_img3_wfss_$basename.sh

