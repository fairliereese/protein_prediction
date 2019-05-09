from pp_utils import pp_utils

gene_names = '/data/users/freese/mortazavi_lab/data/190313_HepG2/190502_HepG2_filtered_protein_pred/blastp/HepG2_filtered_gene_IDS_mini.txt'
p_ref = '/data/users/freese/mortazavi_lab/ref/gencode.v24/gencode.v24.pc_translations.fasta'
odir = '/data/users/freese/mortazavi_lab/data/190313_HepG2/190502_HepG2_filtered_protein_pred/'
prefix = 'HepG2_filtered'
pepfile = '/data/users/freese/mortazavi_lab/data/190313_HepG2/190502_HepG2_filtered_protein_pred/transdecoder/HepG2_filtered.fasta.transdecoder.pep'

b_tbl = pp_utils().run_blastp(gene_names, p_ref, odir, prefix, pepfile)
print(b_tbl)