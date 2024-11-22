#!/usr/bin/env python

#Converts grep data to fasta alignment format. 

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from itertools import islice

counter = 0
file = open(sys.argv[1], 'r')
outfile1 = open(sys.argv[2], 'w+')
outfile2 = open(sys.argv[3], 'w+')

print('output 1: (1) Marginal reconstruction of ancestral sequences \n' +
      '(eqn. 4 in Yang et al. 1995 Genetics 141:1641-1650). \n')

print('output 2: (2) Joint reconstruction of ancestral sequences\n' + 
 '(eqn. 2 in Yang et al. 1995 Genetics 141:1641-1650), \n' +
 'using the algorithm of Pupko et al. (2000 Mol Biol Evol 17:890-896), \n' +
 'modified to generate sub-optimal reconstructions.')

flag = True
for line in file:
	counter += 1
	if 'List of extant and reconstructed sequences' in line.strip() and flag:
		flag = False
		# outfile1.write(line)
		next(file)
		line = file.readline()
		outfile1.write(line)
		no_sequences = int(line.strip().split()[0])
		next(file)
		lines_gen = islice(file, no_sequences)
		for line in lines_gen:
			outfile1.write(line)
	if 'List of extant and reconstructed sequences' in line.strip() and not flag:
		flag = False
		# outfile2.write(line)
		next(file)
		line = file.readline()
		outfile2.write(line)
		no_sequences = int(line.strip().split()[0])
		next(file)
		lines_gen = islice(file, no_sequences)
		for line in lines_gen:
			outfile2.write(line)
file.close()
outfile1.close()
outfile2.close()