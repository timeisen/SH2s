import sys
from Bio import SeqIO

with open(sys.argv[1], 'r') as f1, open(sys.argv[2], 'w+') as f2:
	for record in SeqIO.parse(f1, 'fasta'):
		modseq = str(record.seq).replace("-","")
		f2.write(">" + record.id + "\n" + modseq + "\n")