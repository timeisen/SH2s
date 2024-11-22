from collections import defaultdict as dd
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import regex
import math
R_POS = 33
id_dict = {}
with open(sys.argv[3], 'r') as f: #pid file
	next(f)
	for line in f:
		id_dict[line.split("\t")[0]] = float(line.split("\t")[1].strip())

seq_counter = 0
with open(sys.argv[1], 'r') as f1, open(sys.argv[2], 'w+') as f2: #main code block that parses the fastq file (.extendedFrags )
	for record in SeqIO.parse(f1, "fasta"):
			SeqIO.write(record, f2, 'fasta')
			if record.id in id_dict and id_dict[record.id] > 0.6:
				assert record.seq[R_POS] == 'R'

				seq_new_K = record.seq[:R_POS] + 'K' + record.seq[(R_POS+1):]
				record_K = SeqRecord(seq_new_K, id=record.id + "/CONTROL_K", description = '')
				SeqIO.write(record_K, f2, 'fasta')

				if record.id in id_dict and id_dict[record.id] > 0.95:
					seq_new_X = record.seq[:R_POS] + '*' + record.seq[(R_POS+1):]
					record_X = SeqRecord(seq_new_X, id=record.id + "/CONTROL_X", description = '')
					SeqIO.write(record_X, f2, 'fasta')

				seq_counter += 1

print("seqs mutated:", seq_counter)