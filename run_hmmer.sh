#!/bin/bash


# input arguments
PEP=$1
OUT=$2
TBL=$3

module load hmmer

hmmscan -o $OUT --noali --tblout $TBL /data/users/freese/mortazavi_lab/ref/Pfam-A.hmm $PEP

