
from Bio import SeqIO

# task 2.4
records = SeqIO.parse("globins.gb", "genbank")
for r in records:
    print("{:15s} {:15s} {}".format(r.id, r.name, r.description))
