#!/bin/sh 

#$ -S /bin/sh
#$ -cwd
#$ -V

# # specifing which node
#$ -q all.q@_NODE_

# #$ -l mem_free=120G
# # number of using cores
#$ -pe openmpi 1

# # name of job
#$ -N _JOB_NAME_


# input asn file
asn_file=_ASN_FILE_

# output directory
output_dir=_OUTPUT_DIR_

# log directory
log_dir=_LOG_DIR_

echo 'output dir: '$output_dir
echo 'log dir:    '$log_dir
echo 'asn_file:   '$asn_file

python pipeline_Image3_wfss.py $asn_file $output_dir

mv _JOB_NAME_.o* $log_dir
mv _JOB_NAME_.e* $log_dir
mv exe_pipe_img3__BASENAME_.sh $log_dir
