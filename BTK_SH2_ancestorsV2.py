#BTK_SH2_ancestors.py
#TJE 2022 03 04
#Version2 just uses some BsaI sites. 

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict as dd
import sys
import random 
import json
from types import SimpleNamespace
import pdb
import math


MAX_SEQ_SUM = 1256 #lib size
NUM_ALT_WT = 20
FULL_SEQ_LENGTH = 300
SEQ_REDUNDANCY = 4
PRIMERS = ['CTAGGGAACCAGGCTTAACG', 'GGAAAACTAAGACAAGGCGC', 'TAACGACGTGCCGAACTTAG', 'AGTGAACTGACCGAATCCTC']

# assert False, "I need to use the PFAM alignments otherwise this doesn't work. Including the BTK numbering."

class mutagenesis_region:
	'''defines a class with attributes specific to mutagenesis regions'''
	def __init__(self, wt_dna_seq, region_id, left_handle, right_handle, begin_num, aligned_seq, phiX):
		self.wt_dna_seq = SeqRecord(Seq(wt_dna_seq), id = region_id)
		self.wt_prot = self.wt_dna_seq.translate()
		self.left_handle = left_handle
		self.right_handle = right_handle
		self.phiX = phiX
		self.begin_num = begin_num
		self.aligned_seq = aligned_seq


codon_dict = dd(list)
with open("/Users/timeisen/Dropbox (Personal)/KuriyanLab/Sequences/TwistOrdering/SH2_Variants_20220315/hsap_codon_freq_T.txt", 'r') as f:
	next(f)
	for line in f:
		aa = line.split()[1]
		# if aa == '*': aa = 'X' #new stop codon convention
		codon_dict[aa].append(line.split()[0])

def renumber_btk(pos, region):
	'''Renumbers the aa position by removing the gapped positions in BTK'''
	number_of_gaps = region.aligned_seq[:pos].count('-')
	return(pos - number_of_gaps)		

def populate_seq(region, fasta_handle, codon_dict, NUM_ALT_WT = NUM_ALT_WT):
	'''Main function for reverse translating domains'''
	aa_record, recode_SEQ, recoded_names = [], [], []

	# add in the WT sequence:
	wt_record = SeqRecord(
		region.left_handle + region.wt_dna_seq.seq + region.right_handle,
		id = region.wt_dna_seq.id,
    	name = "",
	    description = "")
	aa_record.append(wt_record)
	#new code as of 2022 06 02
	unknown_aa, gap_counter = 0, 0
	with open(fasta_handle, 'r') as f:
		for record in SeqIO.parse(f, "fasta"): #parsing each seq in alignment file.
			# if 'BTK_HUMAN' record.id
			for i in range(SEQ_REDUNDANCY): #do this n times
				reverse_translate = ''
				for pos, aa in enumerate(str(record.seq).upper()): #make sure that I can iterate over this, or it may need to be converted to a string
					if aa == '-': 
						gap_counter += 1
						continue
					elif aa == region.aligned_seq[pos]: #matches the BTK aa 
						pos_alt = renumber_btk(pos, region)
						reverse_translate += region.wt_dna_seq[(pos_alt*3):(pos_alt*3 + 3)] #the wt dna doesn't include gaps.
					else:
						reverse_translate += random.choice(codon_dict[aa]) #just choose a codon to include at random. 
				NewSeqRecord = SeqRecord(
					Seq(region.left_handle + str(reverse_translate.seq) + region.right_handle),
					id=record.id + "/" + str(i),
    				name="",
	    			description="")
				if 'GGTCTC' in reverse_translate.seq: continue
				if Seq('GGTCTC').reverse_complement() in reverse_translate.seq: continue
				if NewSeqRecord.seq in [record.seq for record in aa_record]: continue
				else: aa_record.append(NewSeqRecord) #
	idx = 0
	wt_fill = MAX_SEQ_SUM - len(aa_record)
	assert wt_fill > NUM_ALT_WT, 'too many sequences'
	while(idx < (wt_fill)):
		recoded_variant, variable_region = recode(region, codon_dict, num_var = 5) #this num_var is the number of codons at which there is variation.
		if recoded_variant not in recode_SEQ: #make sure I don't pick the same ones
			NewWTRecord = SeqRecord(
				Seq(recoded_variant),
				id = region.wt_dna_seq.id + '_' + str(idx),
				name = "",
				description = ""
				)
			if 'GGTCTC' in variable_region: continue
			if str(Seq('GGTCTC').reverse_complement()) in variable_region: continue
			if NewWTRecord.seq in [record.seq for record in aa_record]: continue
			aa_record.append(NewWTRecord)
			idx += 1
	#change asteriks to X for writing


	return(aa_record)

def recode(region, codon_dict, num_var, all = False):
	'''come up with synonomous variants of a particular aa seq'''
	randposaa = random.sample(range(len(region.wt_prot)), k = num_var)
	if all: randposaa = [i for i in range(len(region.wt_prot))]
	newseq = region.wt_dna_seq.seq
	for idx in randposaa:
		wt_codon = region.wt_dna_seq.seq[(idx*3):(idx*3+3)]
		choice_list = [i for i in codon_dict[region.wt_prot[idx]] if i != wt_codon]
		if len(choice_list) < 1: alt_codon = wt_codon
		else: alt_codon = random.choice(choice_list)
		newseq = str(newseq[:(idx * 3)]) + alt_codon + str(newseq[((idx + 1) * 3):])
	return((str(region.left_handle + newseq + region.right_handle), newseq))

def SeqFiller(record, region):
	'''pad sequence'''
	diff = FULL_SEQ_LENGTH - len(record.seq)
	if len(record.seq) >= 280:
		Full_Primer = PRIMERS[0]
		left_pad = Full_Primer[(20 - diff):]
	else:
		Full_Primer = PRIMERS[1] + region.phiX 
		left_pad = Full_Primer[:diff]

	new_seq = left_pad + str(record.seq)
	NewSeqRecord = SeqRecord(
		Seq(new_seq),
		id = record.id,
		name = "",
		description = "")
	if len(NewSeqRecord.seq) != FULL_SEQ_LENGTH:
		print(record.id)
		print(record.seq)
		print(diff)
		assert False, "Sequence length issue"
	return(NewSeqRecord)

def create_seqs(region, fasta_handle, filename):
	'''write files'''
	final_record_list = populate_seq(region, fasta_handle, codon_dict)

	SeqCounter = 0
	with open(filename, 'w+') as f:
		for record in final_record_list:
			if len(record.seq) != FULL_SEQ_LENGTH: record = SeqFiller(record, region)
			SeqIO.write(record, f, 'fasta')
			SeqCounter += 1
	print("Sequences written: ", SeqCounter)

#parse the json file.
def read_json(json_file):
	with open(json_file, 'r') as f:
		all_regions = json.load(f)
		all_regions_class = {}
		for region_name, region in all_regions.items():
			all_regions_class[region_name] = mutagenesis_region(wt_dna_seq = region['wt_dna_seq'],\
				region_id = region['region_id'],\
				left_handle  = region['left_handle'],\
				right_handle  = region['right_handle'],\
				begin_num = region['begin_num'],\
				phiX = region['phiX'],\
				aligned_seq = region['aligned_seq'])
	return(all_regions_class)

all_regions_class = read_json(sys.argv[2]) #json file
create_seqs(all_regions_class['btk_sh2'], sys.argv[1], sys.argv[3]) #fasta handle
