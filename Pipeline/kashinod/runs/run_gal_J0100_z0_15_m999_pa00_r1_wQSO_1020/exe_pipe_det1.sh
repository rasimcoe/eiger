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

# uncal.fits file name
fil=_FILE_NAME_

# output directory
output_dir=_OUTPUT_DIR_

log_dir=_LOG_DIR_

echo 'input file: '$fil
echo 'output dir: '$output_dir 
python pipeline_Detector1.py $fil $output_dir


mv _JOB_NAME_.o* $log_dir
mv _JOB_NAME_.e* $log_dir
mv exe_pipe_det1__BASENAME_.sh $log_dir

echo 'exe_pipe_det1_[job].sh finished.'



