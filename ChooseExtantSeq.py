from collections import defaultdict as dd
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import regex
import math

def pairwise(seq1, seq2):
	id_count = 0
	nid_count = 0
	for pos, aa in enumerate(seq1):
		if aa != '-' and seq2[pos] == aa:
			id_count += 1
		elif aa != '-' and seq2[pos] != aa:
			nid_count += 1
	return(id_count / (id_count + nid_count))

query_dict = {}

with open(sys.argv[1], 'r') as f: #main code block that parses the fasta file for query
	for record in SeqIO.parse(f, "fasta"):
		query_dict[str(record.id)] = str(record.seq)

match_dict = {seqID:[0.0, 'x', 'y'] for seqID, val in query_dict.items()}

with open(sys.argv[2], 'r') as f:
	for record in SeqIO.parse(f, "fasta"):
		for seqID, sequence in query_dict.items():
			id_val = pairwise(sequence, str(record.seq))
			if match_dict[seqID][0] < id_val and str(record.id) not in query_dict.items(): 
				match_dict[seqID][0] = id_val
				match_dict[seqID][1] = str(record.id)
				match_dict[seqID][2] = str(record.seq)


with open(sys.argv[3], 'w+') as f:
	f.write("query\tmatch\tpercent_id\n")
	for query, match_tup in match_dict.items():
		f.write(query + "\t" + match_tup[1] + "\t" + str(match_tup[0]) + "\n")

