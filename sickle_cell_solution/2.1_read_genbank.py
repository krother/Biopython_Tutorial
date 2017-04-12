
from Bio import SeqIO

# task 2.1
records = SeqIO.parse("sickle.gb", "genbank")
print(records)

for rec in records:
    print(rec)
