from collections import defaultdict as dd
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import regex
import math

remove_seqs = []
with open(sys.argv[1], 'r') as f:
	for line in f:
		remove_seqs.append(line.strip())

with open(sys.argv[2], 'r') as f1, open(sys.argv[3], 'w+') as f2: #main code block that parses the fastq file (.extendedFrags )
	for record in SeqIO.parse(f1, "fasta"):
			if str(record.id) in remove_seqs: continue
			else: SeqIO.write(record, f2, 'fasta')