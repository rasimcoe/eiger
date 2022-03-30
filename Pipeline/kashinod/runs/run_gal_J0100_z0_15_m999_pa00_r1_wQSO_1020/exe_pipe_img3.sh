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

# Source catalog arguments (under construction)
# ==> currently, deblending can work only on axion
kernel_fwhm=_KERNEL_FWHM_ # default 2.0
snr_threshold=_SNR_THRESHOLD_ # default 3
npixels=_NPIXELS_ # default 5
deblend=_DEBLEND_ # default False (0)
#echo 'source_catalog arguments::::'
#echo '  kernel_fwhm:      '$kernel_fwhm
#echo '  snr_threshold:    '$snr_threshold
#echo '  npixels:          '$npixels
#echo '  deblend(0=F/1=T): '$deblend

python pipeline_Image3.py $asn_file $output_dir #$kernel_fwhm $snr_threshold $npixels $deblend

mv _JOB_NAME_.o* $log_dir
mv _JOB_NAME_.e* $log_dir
mv exe_pipe_img3__BASENAME_.sh $log_dir

echo 'exe_pipe_img3_[job].sh finished.'

