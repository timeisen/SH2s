#BTK_SH2_ancestors.py
#TJE 2022 03 04

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict as dd
import sys
import random 
import json
from types import SimpleNamespace
import pdb

MAX_SEQ_SUM = 2000 #lib size
NUM_ALT_WT = 40
FULL_SEQ_LENGTH = 300

# assert False, "I need to use the PFAM alignments otherwise this doesn't work. Including the BTK numbering."

class mutagenesis_region:
	'''defines a class with attributes specific to mutagenesis regions'''
	def __init__(self, wt_dna_seq, region_id, left_handle, filler_left, right_handle, filler_right, begin_num, aligned_seq):
		self.wt_dna_seq = SeqRecord(Seq(wt_dna_seq), id = region_id)
		self.wt_prot = self.wt_dna_seq.translate()
		self.left_handle = left_handle
		self.filler_left = filler_left
		self.right_handle = right_handle
		self.filler_right = filler_right
		self.begin_num = begin_num
		self.aligned_seq = aligned_seq


top_codon_dict = dd(list)
with open("/Users/timeisen/Dropbox (Personal)/KuriyanLab/Sequences/TwistOrdering/SH2_Variants_20220315/top_2_codons_hsap_by_usage.txt", 'r') as f:
	next(f)
	for line in f:
		aa = line.split()[1]
		# if aa == '*': aa = 'X' #new stop codon convention
		top_codon_dict[aa].append(line.split()[0])

def renumber_btk(pos, region):
	'''Renumbers the aa position by removing the gapped positions in BTK'''
	number_of_gaps = region.aligned_seq[:pos].count('-')
	return(pos - number_of_gaps)		

def fill_seq(newseq, record, region, full_length = FULL_SEQ_LENGTH):
	'''use flanking regions to bring the length of each seq up to full_length, including handles'''
	distance = full_length - len(newseq) - len(region.left_handle) - len(region.right_handle)
	if distance < 0:
		left_flank = region.left_handle[int(-distance / 2):]
		right_flank = region.right_handle[:int((distance - distance % 2) / 2)] #more on the right side, arbitrarily
		print('trimming handles for:', record.id, distance)
		#THIS IS A PROBLEM STILL, I NEED TO REMOVE THESE SEQS
	else: 
		left_flank = region.filler_left[int(-(distance / 2)):] + region.left_handle
		right_flank = region.right_handle + region.filler_right[:(int((distance / 2)) + distance % 2)] #more sequence on the right side, arbitrarily
	return(left_flank + newseq + right_flank)

def populate_seq(region, fasta_handle, top_codon_dict, NUM_ALT_WT = NUM_ALT_WT):
	'''Main function for reverse translating domains'''
	aa_SEQ, names_SEQ, recode_SEQ, recoded_names = [], [], [], []

	#add in the WT sequence:
	aa_SEQ.append(str(region.left_handle + region.wt_dna_seq.seq + region.right_handle))
	names_SEQ.append(region.wt_dna_seq.id)

	#new code as of 2022 06 02
	unknown_aa, gap_counter = 0, 0
	with open(fasta_handle, 'r') as f:
		for record in SeqIO.parse(f, "fasta"): #parsing each seq in alignment file.
			for i in range(4): #do this n times
				reverse_translate = ''
				for pos, aa in enumerate(str(record.seq)): #make sure that I can iterate over this, or it may need to be converted to a string
					if aa == '-': 
						gap_counter += 1
						continue
					elif aa == region.aligned_seq[pos]: #matches the BTK aa 
						pos_alt = renumber_btk(pos, region)
						reverse_translate += region.wt_dna_seq[(pos_alt*3):(pos_alt*3 + 3)] #the wt dna doesn't include gaps.
					elif aa == 'X':
						pos_alt = renumber_btk(pos, region)
						reverse_translate += region.wt_dna_seq[(pos_alt*3):(pos_alt*3 + 3)]
						unknown_aa += 1
					else:
						reverse_translate += random.choice(top_codon_dict[aa]) #just choose a codon to include at random. 

				newseq_with_flanking = fill_seq(str(reverse_translate.seq), record, region)
				aa_SEQ.append(newseq_with_flanking) #append sequence and name
				# aa_SEQ.append(str(reverse_translate.seq)) #append sequence and name
				names_SEQ.append(record.id + "/" + str(i))
	print('unknown amino acids in alignment:', unknown_aa / 4)
	idx = 0
	while(idx < NUM_ALT_WT):
		recoded_variant = recode(region, top_codon_dict, num_var = 5) #this num_var is the number of codons at which there is variation.
		if recoded_variant not in recode_SEQ: #make sure I don't pick the same ones
			recode_SEQ.append(recoded_variant)
			recoded_names.append(region.wt_dna_seq.id + '_' + str(idx))
			idx += 1
	[aa_SEQ.append(i) for i in recode_SEQ]
	[names_SEQ.append(i) for i in recoded_names]

	#change asteriks to X for writing
	names_SEQ = [val.replace("*", "X") for val in names_SEQ]


	return(names_SEQ, aa_SEQ)

#write file.
def file_writer(file, aa_SEQ, names_SEQ, MAX_SEQ_SUM):
	'''write files'''
	SeqCounter = 0
	with open(file, 'w+') as f:
		for idx, dna_seq in enumerate(aa_SEQ):
			if SeqCounter == MAX_SEQ_SUM: break
			f.write(">TJE_" + names_SEQ[idx] + "\n")
			f.write(dna_seq + "\n")
			SeqCounter += 1
	print("Sequences written: ", SeqCounter)

def recode(region, top_codon_dict, num_var, all = False):
	'''come up with synonomous variants of a particular aa seq'''
	randposaa = random.sample(range(len(region.wt_prot)), k = num_var)
	if all: randposaa = [i for i in range(len(region.wt_prot))]
	newseq = region.wt_dna_seq.seq
	for idx in randposaa:
		wt_codon = region.wt_dna_seq.seq[(idx*3):(idx*3+3)]
		choice_list = [i for i in top_codon_dict[region.wt_prot[idx]] if i != wt_codon]
		if len(choice_list) < 1: alt_codon = wt_codon
		else: alt_codon = random.choice(choice_list)
		newseq = str(newseq[:(idx * 3)]) + alt_codon + str(newseq[((idx + 1) * 3):])
	return(str(region.left_handle + newseq + region.right_handle))



def create_seqs(regionA, fasta_handle, filename):
	names_SEQ_regionA, aa_SEQ_regionA = populate_seq(regionA, fasta_handle, top_codon_dict)

	file_writer(filename, aa_SEQ_regionA, names_SEQ_regionA, MAX_SEQ_SUM)

#parse the json file.
def read_json(json_file):
	with open(json_file, 'r') as f:
		all_regions = json.load(f)
		all_regions_class = {}
		for region_name, region in all_regions.items():
			all_regions_class[region_name] = mutagenesis_region(wt_dna_seq = region['wt_dna_seq'],\
				region_id = region['region_id'],\
				left_handle  = region['left_handle'],\
				filler_left  = region['filler_left'],\
				right_handle  = region['right_handle'],\
				filler_right  = region['filler_right'],\
				begin_num = region['begin_num'],\
				aligned_seq = region['aligned_seq'])
	return(all_regions_class)

all_regions_class = read_json(sys.argv[2]) #json file
create_seqs(all_regions_class['btk_sh2'], sys.argv[1], sys.argv[3]) #fasta handle
