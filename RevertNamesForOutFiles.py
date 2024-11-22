import sys

name_dict = {}
with open(sys.argv[1], 'r') as f:
	next(f) #header
	for line in f:
		original_name = line.split("\t")[0]
		new_name = line.split("\t")[1].strip()
		name_dict[new_name] = original_name

with open(sys.argv[2], 'r') as f1, open(sys.argv[3], 'w+') as f2:
	for line in f1:
		for new_name, original_name in name_dict.items():
			name = line.split()[0]
			if name in name_dict:
				new_name = name_dict[line]
			else: new_name = name
			
		line = line.split(2:)
		line = line.replace(new_name + " ", original_name + " ") #spaces for full seq match
		f2.write(line)