#!/bin/sh

#$ -S /bin/sh
#$ -cwd
#$ -V

# # specifing which node
#$ -q all.q@_NODE_

# # number of using cores
#$ -pe openmpi 1

# # name of job
#$ -N _JOB_NAME_



# input directory
input_dir=_INPUT_DIR_

# output directory
output_dir=_OUTPUT_DIR_

# rate.fits file name
fil=_RATE_FILE_NAME_

log_dir=_LOG_DIR_


echo 'input file: '$fil
echo 'output dir: '$output_dir
python pipeline_Image2.py $fil $output_dir



mv _JOB_NAME_.o* $log_dir
mv _JOB_NAME_.e* $log_dir
mv exe_pipe_img2__BASENAME_.sh $log_dir

echo 'exe_pipe_img2_[job].sh finished.'




