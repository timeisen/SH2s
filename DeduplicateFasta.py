from Bio import SeqIO
import sys
#this deduplicates a fasta file.
#note that I've updated it to retain the SeqID that is for the human sequence. 

record_dict = {}
seq_set = set()
with open(sys.argv[1], "r") as f:
	for record in SeqIO.parse(f, "fasta"):
		if "HUMAN" in record.id:
			record_dict[record.seq] = record
			seq_set.add(record.seq)
		elif record.seq not in seq_set:
			seq_set.add(record.seq)
			record_dict[record.seq] = record

with open(sys.argv[2], "w+") as f:
	for record in record_dict.values():
		SeqIO.write(record, f, "fasta")