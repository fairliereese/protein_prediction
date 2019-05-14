class pp_utils:
	def __init__(self):
		1

	# get basename from a file and path string
	def get_basename(self, filepath):
		import os
		return os.path.basename(os.path.splitext(filepath)[0])

	# get and format output directory
	def format_odir(self, odir):
		import os
		cwd = os.getcwd()

		# if first character is not /, use cwd to make this an absolute path
		if odir[0] != '/' and odir[0] != '~':
		    odir = cwd+odir
		if odir[-1] != '/':
		    odir += '/'
		return odir

	# make a dated output directory for the files used for the tracks
	def make_dated_folder(self, odir, prefix):
		import os
		import datetime

		odir = odir+'/'
		date = datetime.datetime.now()
		date = date.strftime('%y%m%d')
		odir = odir+date+'_'+prefix+'_protein_pred/'

		if not os.path.isdir(odir):
			print('Making new directory '+odir)
			os.makedirs(odir)

		return odir

	# 
	def make_folder(self, odir, dirname):
		import os
		odir = odir+dirname

		if not os.path.isdir(odir):
			print('Making new directory '+odir)
			os.makedirs(odir)

		return odir+'/'

	# get value associated with keyword in the 9th column of gtf
	def get_field_value(self, key, fields):
		import os
		if key not in fields:
			return None
		else:
			return fields.split(key+' "')[1].split()[0].replace('";','')

	#
	def filter_gene_list(self, pepfile, gene_list, odir, prefix):
		import os
		import subprocess

		oname = self.get_pepfile(odir, prefix)

		genefile = open(gene_list, 'r')
		for line in genefile:
			line = line.replace('\n', '')
			genes = line.split(',')
		genefile.close()

		cmd = "grep -A 1 '{}' {} | grep -v '\-\-' > {}".format(\
			'\\|'.join(genes)[:-2], pepfile, oname)
		process = subprocess.call(cmd, shell=True)

		return oname

	#
	def filter_coding_novel_gtf(self, gtffile, odir, prefix):
		import os
		ofile = odir+prefix+'_coding_novel.gtf'

		# print('Making new file '+ofile+'...')

		ofile = open(ofile, 'w')
		infile = open(gtffile, 'r')

		# print(infile.name)
		# print(ofile.name)

		# store gene entry and whether it's been written to ofile
		gene_line = ''
		gene_written = False

		# whether current transcript is protein-coding (or novel)
		t_coding_novel = False

		for line in infile:
			line = line.replace('\n', '')
			temp = line.split('\t')
			fields = temp[-1]

			# gene vs. transcript vs. exon
			if temp[2] == 'transcript':

				# novel transcript
				tstatus = self.get_field_value('transcript_status', fields)
				if tstatus == 'NOVEL': 
					t_coding_novel = True

				else:
					ttype = fields.split('transcript_type "')[1].split()[0]
					ttype = ttype.split('";')[0]

					# protein-coding transcript 
					if ttype == 'protein_coding':
						t_coding_novel = True

					# non-coding transcript
					else: 
						t_coding_novel = False

				# write gene line for transcript if not already written 
				if t_coding_novel == True and gene_written == False:
					ofile.write(gene_line+'\n')
				if t_coding_novel == True:
					ofile.write(line+'\n')

			# store gene line
			elif temp[2] == 'gene':
				gene_line = line
				gene_written = False

			# write exon lines if necessary 
			elif temp[2] == 'exon':
				if t_coding_novel == True:
					ofile.write(line+'\n')

		oname = ofile.name
		ofile.close()
		infile.close()

		return oname

	# generate the gid/tid map
	def gen_tid_gid_map(self, gtffile, prefix, odir):
		import os

		outfile = odir+prefix+'_tid_gid_map.tsv'
		outfile = open(outfile, 'w')

		with open(gtffile, 'r') as gtffile:
			for line in gtffile: 
				temp = line.split('\t')
				if temp[2] == 'transcript':
					fields = temp[8]
					tid = fields.split('transcript_id "')[1].split('";')[0]
					gid = fields.split('gene_id "')[1].split('";')[0]
					outfile.write('\t'.join([gid, tid])+'\n')
		oname = outfile.name	
		outfile.close()
		return oname

	def get_tid_gid_map(self, odir, prefix):
		import os
		prefix = prefix.replace('_gene_list', '')
		f = odir+prefix+'_tid_gid_map.tsv'
		return f

	def run_transdecoder(self, gtffile, fafile, odir, prefix):
		import os
		import subprocess

		bname = self.get_basename(gtffile)

		# make output directory for transdecoder run
		t_odir = self.make_folder(odir, 'transdecoder')
		mapfile = odir+prefix+'_tid_gid_map.tsv'

		# # run transdecoder qsub
		# cmd = "qsub "+\
		# 	   "-v GTF="+gtffile+" "+\
		# 	   "-v BNAME="+bname+" "+\
		# 	   "-v REF="+fafile+" "+\
		# 	   "-v OPATH="+t_odir+" "+\
		# 	   "-v MAP="+map_odir+bname+'_tid_gid_map.tsv'+" "+\
		# 	   "run_transdecoder.sh"

		# run transdecoder
		cmd = 'bash run_transdecoder.sh {} {} {} {} {}'.format(\
			gtffile,\
			prefix,\
			fafile,\
			t_odir,\
			mapfile)
		
		print(cmd)
		process = subprocess.call(cmd, shell=True)

		oname = t_odir+prefix+'.fasta.transdecoder.pep'
		return oname

	def get_pepfile(self, odir, prefix):
		f = odir+'transdecoder/'+prefix+'.fasta.transdecoder.pep'
		return f

	def run_hmmer(self, pepfile, odir, prefix):
		import os
		import subprocess

		bname = self.get_basename(pepfile)

		# make output directory for hmmer run
		odir = self.make_folder(odir, 'hmmer')
		outfile = odir+prefix+'.out'
		outtbl = odir+prefix+'.tab'

		# # run hmmer qsub
		# cmd = "qsub "+\
		# 	"-v PEP="+pepfile+" "+\
		# 	"-v OUT="+bname+" "+\
		# 	"-v TBL="+fafile+" "+\
		# 	"qsub_hmmer.sh"

		# run hmmer 
		cmd = 'bash run_hmmer.sh {} {} {}'.format(\
			pepfile,\
			outfile,\
			outtbl)

		print(cmd)
		process = subprocess.call(cmd, shell=True)

		return outtbl

	# convert to two-line fasta
	def fasta_two_lines(self, pepfile):
		import os

		# write new file
		ofile = open(pepfile+'_temp', 'w')
		with open(pepfile,'r') as f: 
		    i = 0
		    seq = []
		    for line in f:
		        if line[0] == '>' and i == 0:
		            ofile.write(line)
		        elif line[0] == '>':
		            ofile.write(''.join(seq)+'\n')
		            ofile.write(line)
		            seq = []
		        else:
		            seq+=line[:-1]  
		        i+=1
		ofile.write(''.join(seq)+'\n')
		ofile.close()

		# overwrite old file with new data
		os.rename(pepfile+'_temp', pepfile)

	# get tab separated line to write to outfile 
	# different depending on if from reference or from TALON
	def reformat_line(self, line, infile):
		import re
		new_line = []

		temp = line[1:-1].split()
		new_line.append(temp[0].split('.p')[0]) # transcript id
		new_line.append(temp[1].split('~~')[0]) # gene id
		if 'ENS' in new_line[0]: new_line.append(False) # novel transcript
		else: new_line.append(True)
		if 'ENS' in new_line[1]: new_line.append(False) # novel gene
		else: new_line.append(True)
		new_line.append(temp[3].split('type:')[1]) # completeness
		new_line.append(temp[4].split('len:')[1]) # length
		new_line.append(temp[5].split('score=')[1]) # score
		new_line.extend(re.search('([0-9]+)-([0-9]+)', temp[6]).groups()) 
		new_line.append(line[-3]) # strand
		new_line.append(str(infile.tell())) # byte location of sequence
	
		return '\t'.join([str(k) for k in new_line])+'\n'

	# format transdecoder output into a tsv
	def reformat_transdecoder(self, pepfile, prefix):
		import os

		self.fasta_two_lines(pepfile)

		oname = self.format_odir(os.path.dirname(pepfile))+prefix+'.tsv'
		outfile = open(oname, 'w')
		infile = open(pepfile, 'r')

		# write header first
		header = ['transcript_id', 'gene_id', \
         'novel_transcript', 'novel_gene',\
         'completeness', 'length', 'score', 'strand', 
         'start', 'stop', 'seq_byte_loc']
		outfile.write('\t'.join(header)+'\n')

		# loop through file
		line = infile.readline()
		while line:
		    new_line = [] 
		    if line[0] == '>':
		        outfile.write(self.reformat_line(line, infile))
		    line = infile.readline()
		outfile.close()
		infile.close()

		return oname

	# format hmmer output into a tsv
	def reformat_hmmer(self, hmmfile, prefix):
		from Bio import SearchIO
		import os

		oname = self.format_odir(os.path.dirname(hmmfile))+prefix+'.tsv'
		outfile = open(oname, 'w')
		infile = open(hmmfile, 'r')

		# write header to file first
		header = ['transcript_id', 'accession', 'bias', 'bitscore', 'description',\
		         'cluster_num', 'domain_exp_num', 'domain_reported_num',\
		         'env_num', 'evalue', 'id', 'overlap_num', 'region_num']
		outfile.write('\t'.join(header)+'\n')
		for record in SearchIO.parse(infile, 'hmmer3-tab'):
		    
		    tid = record.id.split('|')[0]
		    tid = tid.split('.p')[0]

		    for hit in record.hits:
		        new_line = []
		        new_line.append(tid)
		        new_line.append(str(hit.accession))
		        new_line.append(str(hit.bias))
		        new_line.append(str(hit.bitscore))
		        new_line.append(str(hit.description))
		        new_line.append(str(hit.cluster_num))
		        new_line.append(str(hit.domain_exp_num))
		        new_line.append(str(hit.domain_reported_num))
		        new_line.append(str(hit.env_num))
		        new_line.append(str(hit.evalue))
		        new_line.append(str(hit.id))
		        new_line.append(str(hit.overlap_num))
		        new_line.append(str(hit.region_num))
		        outfile.write('\t'.join(new_line)+'\n')
		infile.close()
		outfile.close()

		return oname

	# get list of gene names associated with the known genes
	def get_gene_names(self, tid_gid_map, odir, prefix, ref_organism):
		import subprocess

		odir = self.make_folder(odir, 'blastp')
		oname = odir+prefix+'_gene_IDS.txt'

		if ref_organism == 'mouse':
			novelty_str = 'ENSMUSG'
		else: novelty_str == 'ENSG'

		cmd = '''grep "{}" {} | cut -f1 | cut -d'.' -f1 > {}'''.format(\
			novelty_str, tid_gid_map, oname)
		print()
		print(cmd)
		print()
		process = subprocess.call(cmd, shell=True)

		return oname

	# run blastp 
	def run_blastp(self, gene_names, p_ref, odir, prefix, pepfile):
		import subprocess 

		odir = odir+'blastp/'
		tab_oname = odir+prefix+'_blast_results.tab'
		cmd = 'bash run_blastp.sh {} {} {} {} {}'.format(\
			gene_names,\
			p_ref,\
			odir,\
			prefix,
			pepfile)
		print()
		print(cmd)
		print()
		subprocess.call(cmd, shell=True)

		return tab_oname

	# reformat blastp output
	def reformat_blastp(self, blastfile, prefix):
		import os

		oname = self.format_odir(os.path.dirname(blastfile))+prefix+'.tsv'
		ofile = open(oname, 'w')
		infile = open(blastfile, 'r')

		#
		header = ['query_tid', 'novel_transcript', 'ref_tid',\
			'ref_gid', 'bitscore',\
			'percent_identity', 'evalue', 'mismatches',\
			'ref_start', 'ref_end', 'query_start', 'query_end',\
			'align_len', 'query_len', 'subject_len', 'num_gaps']
		ofile.write('\t'.join(header)+'\n')

		for line in infile:
			new_line = [] 
			line = line.split('\t')

			acc = line[0]
			# query transcript id
			tid = acc.split('|')[0]
			tid = tid.split('.p')[0]
			new_line.append(tid) 
			# transcript novelty
			if 'ENS' in tid: new_line.append('False')
			else: new_line.append('True')
			# subject transcript id
			new_line.append(line[1].split('|')[1])
			# subject gene id
			new_line.append(line[1].split('|')[2])

			# append the rest of the alignment information
			new_line+=line[2:]
			ofile.write('\t'.join([str(x) for x in new_line]))

		ofile.close()
		infile.close()

		return oname

	#
	def gen_report(self, t_tsv, h_tsv, b_tsv, odir, prefix):
		import pandas as pd

		odir = self.make_folder(odir, 'figures')

		t_df = pd.read_csv(t_tsv, '\t')
		h_df = pd.read_csv(h_tsv, '\t')
		b_df = pd.read_csv(b_tsv, '\t')

		# how many transcripts match up with the reference known transcript
		ofile = odir+prefix+'_known_tid_match.csv'
		with open(ofile, 'w') as ofile:
			known = b_df.loc[b_df.novel_transcript == False]
			n_known = len(known.query_tid.unique())
			n_match = len(known.loc[known_blast.query_tid == known.ref_tid])
			n_nomatch = abs(n_known-n_match)

			ofile.write(','.join(['Novelty', 'Count', 'Match'])+'\n')
			ofile.write(','.join(['Known', str(n_match), 'Matches'])+'\n')
			ofile.write(','.join(['Known', str(n_nomatch), 'Does not match'])+'\n')

		# how similar is each new protein isoform to existing isoforms (blast)
		ofile = odir+prefix+'_novel_blast_match_pctg.csv'
		with open(ofile, 'w') as ofile:
			n_blast = b_df.loc[b_df.novel_transcript == True]\
						.groupby('query_tid').first().reset_index()
			n_blast.to_csv(ofile, index=False)

		# domains in each transcript
		ofile = odir+prefix+'_domains_per_transcript.csv'
		with open(ofile, 'w') as ofile:
			temp_df = h_df.groupby('transcript_id').size().to_frame()
			temp_df['Novelty'] = temp_df.apply(lambda x: \
									 'Known' if 'ENS' in x.index \
									  else 'Novel', axis=1)
			temp_df = temp_df.rename({0: 'Count'}, axis=1)
			temp_df = pd.merge(temp_df, t_df, on='transcript_id')
			temp_df = temp_df.rename({'completeness': 'Completeness'}, axis=1)
			temp_df.to_csv(ofile, ',', index=False)













