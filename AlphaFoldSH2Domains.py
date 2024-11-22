#load library
from pymol import cmd, stored
import sys
import os
#cd /Users/timeisen/Dropbox (Personal)/KuriyanLab/Sequences/TwistOrdering/SH2_Variants_20220315
filenames = '/Users/timeisen/Dropbox (Personal)/KuriyanLab/Alphafold/PFAM_SH2_filenames_UP000005640_9606_HUMAN_v2_filenames_cleaned.txt'
Ids = '/Users/timeisen/Dropbox (Personal)/KuriyanLab/Sequences/AncestralSequenceReconstruction/SH2_Domain_Reconstruction/sh2_reconstruction_2/PF00017_HUMAN_IDs_filtered.txt'
Uniprot_match = '/Users/timeisen/Dropbox (Personal)/KuriyanLab/Sequences/AncestralSequenceReconstruction/SH2_Domain_Reconstruction/sh2_reconstruction_2/PF00017_HUMAN_SH2_proteins_Uniprot_IDs_with_names.txt'
id_dict, Uniprot_match_dict = {}, {}
loaded_files = set()
loaded_files_human_names = set()
with open(Ids, 'r') as f:
	next(f)
	for line in f:
		id_dict[line.split("\t")[0].split("/")[0]] = line.split("\t")[0].split("/")[1]

with open(Uniprot_match, 'r') as f:
	for line in f:
		#pymol name = [number, human_name]
		Uniprot_match_dict["AF-" + line.split("\t")[1].strip() + "-F1-model_v2"] = [id_dict[line.split("\t")[0]], line.split("\t")[0]]

with open(filenames, 'r') as f:
	for filename in f:
		loaded_files.add(filename.strip()[:-4])
		cmd.load('/Users/timeisen/Dropbox (Personal)/KuriyanLab/Alphafold/Human_SH2_domain_proteins/UP000005640_9606_HUMAN_v2/' + filename.strip())
		
for pymol_name, [resi_region, human_name] in Uniprot_match_dict.items():
	if pymol_name not in loaded_files: continue #not all of the uniprot proteins are in the alphafold database
	# print(human_name, "resi" + " " + resi_region + " and " + pymol_name)
	loaded_files_human_names.add(human_name)
	cmd.create(human_name, "resi" + " " + resi_region + " and " + pymol_name)
	cmd.delete(pymol_name)

for human_name in loaded_files_human_names:
	if human_name == 'BTK_HUMAN': continue #doesn't work with btk matching btk
	cmd.super(human_name, "BTK_HUMAN")