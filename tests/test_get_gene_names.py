from pp_utils import pp_utils

tid_gid_map = '/Users/fairliereese/mortazavi_lab/data/190313_HepG2/HepG2_filtered_talon_coding_novel_tid_gid_map.tsv'
odir = '/Users/fairliereese/mortazavi_lab/data/190313_HepG2/190502_HepG2_filtered_protein_pred/'
prefix = 'HepG2_filtered'
gene_names = pp_utils().get_gene_names(tid_gid_map, odir, prefix)
print(gene_names)