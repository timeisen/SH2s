from collections import defaultdict as dd
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import regex
import math
from Bio.Restriction import *

counter = 0
all_seqs = 0
with open(sys.argv[1], 'r') as f1: #main code block that parses the fastq file (.extendedFrags )
	for record in SeqIO.parse(f1, "fasta"):
		all_seqs += 1
		if len(XmaI.search(record.seq)) > 0: 
			print('XmaI', record.id)
			counter += 1
		if len(BamHI.search(record.seq)) > 0: 
			print('BamHI', record.id)
			counter += 1

print(counter)
print(all_seqs)