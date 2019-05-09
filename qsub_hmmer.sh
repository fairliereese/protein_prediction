#!/bin/bash
#$ -N transdecoder
#$ -pe openmp 8-64
#$ -cwd
#$ -q free64,bio,som
#$ -e ~/mortazavi_lab/qsub_output
#$ -o ~/mortazavi_lab/qsub_output
#$ -ckpt restart

# input arguments
# PEP
# OUT
# TBL

module load hmmer

# hmmscan -o ${OPATH}${PREFIX}hmmer.out --noali --tblout ${OPATH}${PREFIX}hmmer_table.out /data/users/freese/mortazavi_lab/ref/Pfam-A.hmm ${DPATH}HepG2_filtered_talon_td.fasta.transdecoder.pep
hmmscan -o $OUT --noali --tblout $TBL /data/users/freese/mortazavi_lab/ref/Pfam-A.hmm $PEP

