
from Bio import SeqIO

# task 2.2
records = SeqIO.parse("sickle.gb", "genbank")
for rec in records:
    print(dir(rec))
    print(rec.id)
    print(rec.name)
    print(rec.description)
