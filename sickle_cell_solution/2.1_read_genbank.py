
from Bio import SeqIO

# task 2.1
records = SeqIO.parse("sickle.gb", "genbank")
print type(records)
print dir(records)
	
