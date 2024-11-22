from collections import defaultdict as dd
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import regex
import math

#arg1 is fasta alignment, arg2 is tab-delim output
with open(sys.argv[1], 'r') as f1, open(sys.argv[2], 'w+') as f2: #main code block that parses the fastq file (.extendedFrags )
	f2.write('id\tpos\taa\n')
	for record in SeqIO.parse(f1, "fasta"):
		for i, aa in enumerate(record.seq):
			newline = record.name + '\t' + str(i) + '\t' + str(aa) + '\n'
			f2.write(newline)