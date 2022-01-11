#!/bin/sh

#$ -S /bin/sh
#$ -cwd
#$ -V

# # specifing which node
#$ -q all.q@messier02

# # number of using cores
#$ -pe openmpi 1

# # name of job
#$ -N gen_yaml

python generate_yaml_files.py

