from collections import defaultdict as dd
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import regex
import math
from Bio.Restriction import *
from Bio.SeqIO import FastaIO


counter = 0
all_seqs = 0
with open(sys.argv[1], 'r') as f1, open(sys.argv[2], 'w+') as f2: #main code block that parses the fastq file (.extendedFrags )
	fasta_out = FastaIO.FastaWriter(f2, wrap=None)
	for record in SeqIO.parse(f1, "fasta"):
		if counter < 2000:
			writeseq = SeqIO.SeqRecord(BsaI.catalyse(record.seq)[0], record.id, description = '')
		else: 
			writeseq = SeqIO.SeqRecord(BsaI.catalyse(record.seq)[1], record.id, description = '')
		counter += 1
		# SeqIO.write(writeseq, f2, 'fasta', width = 1000)
		fasta_out.write_record(writeseq)

