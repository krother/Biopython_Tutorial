
from Bio import SeqIO

# task 2.2
records = SeqIO.parse("sickle.gb", "genbank")
for r in records:
    print r.id
    print r.name
    print r.description
