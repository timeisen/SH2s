from collections import defaultdict as dd
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import regex
import math

BTK = sys.argv[3] #BTK_seq

def pairwise(seq1, seq2):
	seq1 = seq1.upper() #ignore case
	seq2 = seq2.upper()
	id_count = 0
	nid_count = 0
	for pos, aa in enumerate(seq1):
		if aa != '-' and seq2[pos] == aa:
			id_count += 1
		elif aa != '-' and seq2[pos] != aa:
			nid_count += 1
	return(id_count / (id_count + nid_count))

id_dict = {}

with open(sys.argv[1], 'r') as f: #main code block that parses the fastq file (.extendedFrags )
	for record in SeqIO.parse(f, "fasta"):
		pairwise_identity = pairwise(BTK, str(record.seq))
		id_dict[record.id] = pairwise_identity

with open(sys.argv[2], 'w+') as f:
	f.write('Seq_ID\tpairwise_identity\n')
	for seq, pairwise_identity in id_dict.items():
		f.write(seq + "\t" + str(pairwise_identity) + "\n")