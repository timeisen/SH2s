from collections import defaultdict as dd
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import regex
import math

BTK = "---------------------------------------------WYS-K--H-M---T--R---------------------SQ--A-E-Q-L-LKQ---------------------EGK--E--G-G--FI-------------V--R-----------D------S--------S--------------K----------A--------G-----------K---------Y-------T----V---S--VFA--KStgd--------------------------------------pqgVIR-H----Y--V----V--C-----S---T---P--QS--------------------------------Q-Y-Y---L------A-----E------K-----HL-------F--S--T--I-P-ELINYH------------------------------------------------------------------------------------------"

def pairwise(seq1, seq2):
	id_count = 0
	for pos, aa in enumerate(seq1):
		if aa != '-' and seq2[pos] == aa:
			id_count += 1
	return(id_count / 82)

id_dict = {}
linspace_pairwise = {}

linspace_pairwise_counters = {}
for i in range(3, 11):
	linspace_pairwise_counters[i / 10] = 50
seq_list = []
fasta_writer_list = []
with open(sys.argv[1], 'r') as f: #main code block that parses the fastq file (.extendedFrags )
	for record in SeqIO.parse(f, "fasta"):
		pairwise_identity = pairwise(BTK, str(record.seq))
		id_dict[record.id] = pairwise_identity
		if (math.floor(pairwise_identity * 10) / 10) < .3: continue
		elif linspace_pairwise_counters[math.floor(pairwise_identity * 10) / 10] > 0 and str(record.seq) not in seq_list:
			linspace_pairwise[record.id] = pairwise_identity
			linspace_pairwise_counters[math.floor(pairwise_identity * 10) / 10] -= 1
			seq_list.append(str(record.seq))
			fasta_writer_list.append(record)

# with open(sys.argv[2], 'w+') as f:
# 	f.write('Seq_ID\tpairwise_identity\n')
# 	for seq, pairwise_identity in id_dict.items():
# 		f.write(seq + "\t" + str(pairwise_identity) + "\n")

with open(sys.argv[2], 'w+') as f:
	f.write('Seq_ID\tpairwise_identity_arrayed\n')
	for seq, pairwise_identity in linspace_pairwise.items():
		f.write(seq + "\t" + str(pairwise_identity) + "\n")
# 

with open(sys.argv[3], 'w+') as f:
	SeqIO.write(fasta_writer_list, f, 'fasta')