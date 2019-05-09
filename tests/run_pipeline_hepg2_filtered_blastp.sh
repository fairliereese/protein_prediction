#!/bin/bash
#$ -N ppp
#$ -pe openmp 8-64
#$ -cwd
#$ -q free64,bio,som
#$ -e ~/mortazavi_lab/qsub_output
#$ -o ~/mortazavi_lab/qsub_output
#$ -ckpt restart

GTF=/data/users/freese/mortazavi_lab/data/190313_HepG2/HepG2_filtered_talon.gtf
PREFIX=HepG2
ODIR=/data/users/freese/mortazavi_lab/data/190313_HepG2/190502_HepG2_filtered_protein_pred/
OPTS='--filter_noncoding --blastp'

python protein_pred_pipeline.py --gtf $GTF --prefix $PREFIX  --odir $ODIR $OPTS