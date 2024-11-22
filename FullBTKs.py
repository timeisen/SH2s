from Bio import SeqIO
from Bio.Seq import Seq
import sys


btk_records = {}
with open('/Users/timeisen/Dropbox (Personal)/KuriyanLab/Sequences/AncestralSequenceReconstruction/SH2_Domain_Reconstruction/sh2_reconstruction_2/btk_handles.txt', 'r') as f:
	for record in SeqIO.parse(f, 'fasta'):
		btk_records[record.id] = record

record_list = []
with open(sys.argv[1], 'r') as f:
	for record in SeqIO.parse(f, 'fasta'):
		record.seq = Seq(str(record.seq).replace('-', '').upper())
		record.seq = btk_records['BTK_left_no_PHTH_209-280'].seq + record.seq + btk_records['BTK_right_363-660'].seq

		record_list.append(record)

with open(sys.argv[2], 'w+') as f:
	for record in record_list:
		SeqIO.write(record, f, 'fasta')
