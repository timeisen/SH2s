#Makes a unique and short name for every fasta entry.

import sys, argparse
from collections import defaultdict as dd
from Bio import SeqIO
import re
from Bio.SeqRecord import SeqRecord

output_table = open(sys.argv[2],'w+')
output_fasta = open(sys.argv[3],'w+')

output_table.write("original_name\tnew_name\n")

counter = 1 #matching 1-base codeml ordering
for seq_record in SeqIO.parse(sys.argv[1], "fasta"):
	original_name = seq_record.name
	new_name = 'Seq_' + str(counter)
	new_record = SeqRecord(
    	seq_record.seq,
    	id = new_name,
    	description="")
	counter += 1
	output_table.write(original_name + '\t' + new_name + '\n')
	SeqIO.write(new_record, output_fasta, "fasta")
output_table.close()
output_fasta.close()
