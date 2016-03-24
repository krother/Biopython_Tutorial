
from Bio import SeqIO

# task 2.4
records = SeqIO.parse("globins.gb", "genbank")
for r in records:
    print r.id, 
    print r.name, 
    print r.description
