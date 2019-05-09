#!/bin/bash


# input arguments
PEP=$1
OUT=$2
TBL=$3

module load hmmer

# hmmscan -o ${OPATH}${PREFIX}hmmer.out --noali --tblout ${OPATH}${PREFIX}hmmer_table.out /data/users/freese/mortazavi_lab/ref/Pfam-A.hmm ${DPATH}HepG2_filtered_talon_td.fasta.transdecoder.pep
hmmscan -o $OUT --noali --tblout $TBL /data/users/freese/mortazavi_lab/ref/Pfam-A.hmm $PEP

