#!/bin/bash
#$ -N ppp
#$ -pe openmp 8-64
#$ -cwd
#$ -q free64,bio,som
#$ -e ~/mortazavi_lab/qsub_output
#$ -o ~/mortazavi_lab/qsub_output
#$ -ckpt restart

python protein_pred_pipeline.py