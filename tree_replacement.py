#python code to replace a tree with different names
import sys

replace_dict = {} 
with open(sys.argv[1], 'r') as f:
	next(f)
	for line in f:
		replace_dict[line.split("\t")[0]] = line.split("\t")[1].strip()


output = open(sys.argv[3], 'w+')
with open(sys.argv[2], 'r') as f:
	for line in f:
		for key, value in replace_dict.items():
			line = line.replace(key, value)
		output.write(line)
output.close