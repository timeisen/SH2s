from collections import defaultdict as dd
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import regex
import math

with open(sys.argv[1], 'r') as f1, open(sys.argv[2], 'w+') as f2: #main code block that parses the fastq file (.extendedFrags )
	for record in SeqIO.parse(f1, "fasta"):
			if len(str(record.seq).replace("-","")) <= int(sys.argv[3]): SeqIO.write(record, f2, 'fasta')
			else: continue