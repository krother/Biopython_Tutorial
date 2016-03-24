
from Bio import SeqIO

# task 2.3
records = SeqIO.parse("sickle.gb", "genbank")
SeqIO.write(records, "sickle.fasta", "fasta")

