#!/bin/bash
#$ -N transdecoder
#$ -pe openmp 8-64
#$ -cwd
#$ -q free64,bio,som
#$ -e ~/mortazavi_lab/qsub_output
#$ -o ~/mortazavi_lab/qsub_output
#$ -ckpt restart

echo ${GTF}
echo ${BNAME}
echo ${REF}
echo ${OPATH}
echo ${MAP}

UTILPATH=/data/apps/user_contributed_software/freese/TransDecoder/5.5.0/util/

# make the tid/gid mapping
python ~/mortazavi_lab/bin/gen_tid_gid_map_gtf.py --g $GTF

export PERL_HASH_SEED=0
cd $OPATH
module load freese/TransDecoder

echo 'Generating fasta sequences from GTF'
${UTILPATH}gtf_genome_to_cdna_fasta.pl $GTF $REF > ${OPATH}${BNAME}.fasta
${UTILPATH}gtf_to_alignment_gff3.pl $GTF > ${OPATH}${BNAME}.gff3

echo 'Running LongOrfs'
TransDecoder.LongOrfs --gene_trans_map $MAP -S -t ${OPATH}${BNAME}.fasta
echo 'Running Predict'
TransDecoder.Predict --single_best_only -t ${OPATH}${BNAME}.fasta

${UTILPATH}cdna_alignment_orf_to_genome_orf.pl ${OPATH}${BNAME}.fasta.transdecoder.gff3 ${OPATH}${BNAME}.gff3 ${OPATH}${BNAME}.fasta > ${OPATH}${BNAME}.fasta.transdecoder.genome.gff3
