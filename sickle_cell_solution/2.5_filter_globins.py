
from Bio import SeqIO

# task 2.5
records = SeqIO.parse("globins.gb", "genbank")
for r in records:
    if not ('vector' in r.description or 'isolate' in r.description):
        print("{:15s} {:15s} {}".format(r.id, r.name, r.description))
