
from Bio import SeqIO

# task 3.2
records = SeqIO.parse("sickle.gb", "genbank")
sickle = records.next()
rna = sickle.seq.transcribe()
print rna
