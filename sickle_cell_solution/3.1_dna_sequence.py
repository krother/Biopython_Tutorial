
from Bio import SeqIO

# task 3.1
records = SeqIO.parse("sickle.gb", "genbank")
sickle = next(records)
print(dir(sickle))

print(sickle.id)
print(sickle.description)
print(sickle.seq)
