from collections import defaultdict as dd
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import regex
import math
from operator import itemgetter

SEQ_REDUNDANCY = int(sys.argv[4])

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

match_dict = {seqID:[] for seqID, val in query_dict.items()}

all_match_keys = set()
[all_match_keys.add(x) for x in query_dict.keys()] #don't find any of those seqs

total_seq_counter = 0
with open(sys.argv[2], 'r') as f:
	for record in SeqIO.parse(f, "fasta"):
		total_seq_counter += 1
		if total_seq_counter % 10000 == 0: 
			print("seqs processed:")
			print(total_seq_counter)
		for seqID, sequence in query_dict.items():
			id_val = pairwise(sequence, str(record.seq))
			if match_dict and str(record.id) not in all_match_keys: 
				match_dict[seqID].append([id_val, record])
				# all_match_keys.add(str(record.id))



new_match_dict = dd(list)
for query, found in match_dict.items():
	counter = 0
	for entry in sorted(found, key=itemgetter(0), reverse = True): #how many seqs are retained?
		if entry[0] >= 0.987 or flag: 
			flag = False
			continue #single point mutants. 
		elif counter < SEQ_REDUNDANCY: 
			new_match_dict[query].append(entry)
			flag = True #every other sequence. 
			counter += 1
		else: break

write_set = set()
with open(sys.argv[3], 'w+') as f:
	f.write("query\tmatch\tpercent_id\n")
	for query, match_list in new_match_dict.items(): #a dictionary of list of lists, Query_ID: [[Pairwise_val, Subj_record], [Pairwise_val, Subj_ID], ...]
		for record_elem in match_list:
			f.write(query + "\t" + record_elem[1].id + "\t" + str(record_elem[0]) + "\n") 
			##I used to do all of this fancy stuff to remove redundancy, but now I just print it and deal with it in unix. I think it's better to have the best alignments in the file anyway.

			# if record_elem[1].id in write_set and SEQ_REDUNDANCY > 0:
			# 	SEQ_REDUNDANCY -= 1
			# elif record_elem[1].id not in write_set:
			# 	f.write(query + "\t" + record_elem[1].id + "\t" + str(record_elem[0]) + "\n")
			# 	write_set.add(record_elem[1].id)
			# else: 
			# 	f.write(query + "\n")