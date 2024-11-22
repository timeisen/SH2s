from collections import defaultdict as dd
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import regex
import math


def pairwise(seq1, seq2):
	id_count = 0
	nid_count = 0
	for pos, aa in enumerate(seq1):
		if aa != '-' and seq2[pos] == aa:
			id_count += 1
		elif aa != '-':
			nid_count += 1
	return(id_count / (id_count + nid_count))

id_set = set()
with open(sys.argv[1], 'r') as f: #main code block that parses the fastq file (.extendedFrags )
	for record in SeqIO.parse(f, "fasta"):
		id_set.add((str(record.seq), str(record.id)))

pairwise_arr = np.empty(np.array([len(id_set),len(id_set)]), dtype = 'float64')

#create matrix and write column names
with open(sys.argv[3], 'w+') as f:
	f.write("col_names\n")
	for i, seq1_tup in enumerate(id_set):
		f.write(seq1_tup[1] + "\n")
		for j, seq2_tup in enumerate(id_set):
			pairwise_arr[i, j] = pairwise(seq1_tup[0], seq2_tup[0])

np.savetxt(sys.argv[2], pairwise_arr, delimiter = '\t') #save the matrix
