
from Bio import SeqIO

# task 3.1
records = SeqIO.parse("sickle.gb", "genbank")
sickle = records.next()
print dir(sickle)

print sickle.id
print sickle.description
print sickle.seq

