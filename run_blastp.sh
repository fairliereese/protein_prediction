#!/bin/bash

# input arguments
GNAMES=$1
PREF=$2
DIR=$3
PREFIX=$4
PEP=$5

module load blast

# remove contents from output files
> ${DIR}blast.log
> ${DIR}db_creation.log

set +x

lines=`cat ${GNAMES}`
for line in $lines;
do
  # grep valid entries from reference and query
  grep -A 1 $line ${PREF} | grep -v -- "^--$" > ${DIR}${line}_ref.fasta
  grep -A 1 $line ${PEP} | grep -v -- "^--$" > ${DIR}${line}_query.fasta

  # make blast database for that gene
  makeblastdb -in ${DIR}${line}_ref.fasta -dbtype "prot" &>> ${DIR}db_creation.log

  echo $line >> ${DIR}blast.log
  
  # blast query sequences against reference sequences
  blastp -query ${DIR}${line}_query.fasta -db ${DIR}${line}_ref.fasta -out ${DIR}${line}.out -num_alignments 5 -outfmt '6 qacc sallseqid bitscore pident evalue mismatch sstart send qstart qend length qlen slen gaps' -max_hsps 1 &>> ${DIR}${PREFIX}_blast.out
 
  # stack everything in same file
  cat ${DIR}${line}.out >> ${DIR}${PREFIX}_blast_results.tab
  
  echo -e '\n' >> ${DIR}${PREFIX}_blast.out

  # clean up
  rm ${DIR}${line}_ref.fasta
  rm ${DIR}${line}_ref.fasta.phr
  rm ${DIR}${line}_ref.fasta.pin
  rm ${DIR}${line}_ref.fasta.psq
  rm ${DIR}${line}_query.fasta
  rm ${DIR}${line}.out
done

