import sys

name_dict = {}
with open(sys.argv[1], 'r') as f:
	next(f) #header
	for line in f:
		original_name = line.split("\t")[0]
		new_name = line.split("\t")[1].strip()
		new_name_alt = new_name.split("_")[-1] + "_" + new_name
		name_dict[new_name_alt] = original_name

with open(sys.argv[2], 'r') as f1, open(sys.argv[3], 'w+') as f2:
	for line in f1:
		for new_name_alt, original_name in name_dict.items():
			line = line.replace(new_name_alt, original_name)
		f2.write(line)