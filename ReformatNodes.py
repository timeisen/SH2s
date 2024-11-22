import sys
#takes an out1 file format and convers the nodes to a fasta file

with open(sys.argv[1], 'r') as f1, open(sys.argv[2], 'w+') as f2:
	for line in f1:
		if line.split()[0] != 'node': continue
		line_id = line.split()[0] + line.split()[1]
		seq = "".join(line.split()[2:])
		f2.write(">" + line_id + "\n" + seq + "\n")