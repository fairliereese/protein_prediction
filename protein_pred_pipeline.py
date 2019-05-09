# protein_pred_pipeline.py

from pp_utils import pp_utils
import os
from optparse import OptionParser

# # HepG2 HPC
# prefix = 'HepG2'
# gtffile = '/data/users/freese/mortazavi_lab/data/190313_HepG2/HepG2_filtered_talon.gtf' # filtered for bioreps
# fastafile = '/data/users/freese/mortazavi_lab/ref/hg38/hg38.fa'
# p_ref = '/data/users/freese/mortazavi_lab/'

# things_to_run = [\
# 				 'filter_noncoding',\
# 				 'transdecoder',\
# 				 'hmmer',\
# 				 'blastp',\
# 				 'plotting',\
# 				 ]

# arg parsing
parser = OptionParser()
parser.add_option('--gtf', '-g', dest='gtffile',
	help='GTF file to predict proteins from', type='string')
parser.add_option('--prefix', dest='prefix', 
	help='prefix for output directory name. do not include path info',
	type='string')
parser.add_option('--odir', dest='odir', 
	help='output directory IF YOU HAVE ALREADY RUN PART of this program.\
	      do not choose this option if you are starting a new run.',
	      type='string', default=None)
parser.add_option('--filter_noncoding', dest='filter_noncoding',
	action='store_true', help='remove known, non-coding transcripts\
		    before running any protein stuff', default=None)
parser.add_option('--transdecoder', dest='transdecoder',
	action='store_true', help='Run transdecoder on GTF', default=None)
parser.add_option('--hmmer', dest='hmmer', action='store_true',
	help='run hmmer on transdecoder output. requires --transdecoder or --odir \
		  with preexisiting transdecoder data', default=None)
parser.add_option('--blastp', dest='blastp', action='store_true',
	help='run blastp on transdecoder output, requires --transdecoder or --odir \
		  with preexisting transdecoder data', default=None)
parser.add_option('--report', dest='report', action='store_true', 
	help='generate tables and figures associated with output. requires \
		  --transdecoder, --hmmer, and --blastp or --odir with \
		  preexisting data from all 3', default=None)
parser.add_option('--gene_list', dest='gene_list', default=None,\
	help='file with comma-separated list of gene ensembl ids to \
		  specifically examine')
parser.add_option('--ref_organism', dest='ref_organism', default='human', 
	help='organism reference to use. default is human. option is mouse.')
(options, _) = parser.parse_args()

gtffile = options.gtffile
prefix = options.prefix
odir_opt = options.odir
filter_noncoding = options.filter_noncoding
transdecoder = options.transdecoder
hmmer = options.hmmer
blastp = options.blastp
report = options.report 
gene_list = options.gene_list
ref_organism = options.ref_organism

# get parent output file that everything's gonna go into 
if not odir_opt:
	odir = pp_utils().make_dated_folder(os.path.dirname(gtffile), prefix)
else: 
	odir = odir_opt

# set reference files based on what organism we're using
# reference files that we'll keep hard coded for now
if ref_organism == 'human':
	fastafile = '/data/users/freese/mortazavi_lab/ref/hg38/hg38.fa'
	p_ref = '/data/users/freese/mortazavi_lab/ref/gencode.v24/gencode.v24.pc_translations.fasta'
elif ref_organism == 'mouse':
	fastafile = '/data/users/freese/mortazavi_lab/ref/mm10/mm10.fa'
	p_ref = '/data/users/freese/mortazavi_lab/ref/gencode.vM21/gencode.vM21.pc_translations.fasta'

# only analyze genes that are in a given list
if gene_list:
	prefix = prefix+'_gene_list'
	gtffile = pp_utils().filter_gene_list(gtffile, gene_list, odir, prefix)


# filter out known non-coding transcripts from talon gtf
if filter_noncoding:
	prefix = prefix+'_filtered'
	gtffile = pp_utils().filter_coding_novel_gtf(gtffile, odir, prefix)

# run transdecoder
if transdecoder:
	tid_gid_map = pp_utils().gen_tid_gid_map(gtffile, prefix, odir) 
	pepfile = pp_utils().run_transdecoder(gtffile, fastafile, odir, prefix) 
	t_tsv = pp_utils().reformat_transdecoder(pepfile, prefix)

# run hmmer
if hmmer:
	if not transdecoder and not odir_opt:
		print('You must run transdecoder or provide a directory \
			   with preexisting transdecoder data to run hmmer')
		quit()
	hmm_tbl = pp_utils().run_hmmer(pepfile, prefix)
	h_tsv = pp_utils().reformat_hmmer(hmm_tbl, prefix)

# run blastp
if blastp:
	if not transdecoder and not odir_opt:
		print('You must run transdecoder or provide a directory \
			   with preexisting transdecoder data to run blastp')
		quit()
	elif not transdecoder and odir_opt:
		tid_gid_map = pp_utils().get_tid_gid_map(odir, prefix)
		print('Using tid gid map file: '+tid_gid_map)
		pepfile = pp_utils().get_pepfile(odir, prefix)
		print('Using pepfile: '+pepfile)
		
	gene_names = pp_utils().get_gene_names(tid_gid_map, odir, prefix)
	b_tbl = pp_utils().run_blastp(gene_names, p_ref, odir, prefix, pepfile)
	b_tsv = pp_utils().reformat_blastp(b_tbl, prefix)

# plot some cool stuff
if report:
	if not transdecoder or not hmmer or not blastp:
		if not odir_opt:
			print('You must run transdecoder, hmmer, and blastp or \
				   provide a directory will all such preexisting data \
				   to make a report.')
			quit()
	1




