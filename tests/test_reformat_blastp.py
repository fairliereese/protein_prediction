from pp_utils import pp_utils

blastfile = '/Users/fairliereese/mortazavi_lab/data/190313_HepG2/190502_HepG2_filtered_protein_pred/blastp/HepG2_filtered_blast_results.tab'
prefix = 'HepG2_filtered'

b_tsv = pp_utils().reformat_blastp(blastfile, prefix)
print(b_tsv)