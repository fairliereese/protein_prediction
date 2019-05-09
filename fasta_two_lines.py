import argparse
import os

# argument parsing
parser = argparse.ArgumentParser(description='Reformats fasta files that have sequences on more than one line')
parser.add_argument('--f', help='fasta file to convert')
args = parser.parse_args()
ffile = args.f

# write new file
ofile = open(ffile+'_temp', 'w')
with open(ffile,'r') as f: 
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
os.rename(ffile+'_temp', ffile)

