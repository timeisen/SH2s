from Bio import SeqIO
import sys

records = SeqIO.parse(sys.argv[1], "stockholm")
count = SeqIO.write(records, sys.argv[2], "fasta")
print("Converted %i records" % count)