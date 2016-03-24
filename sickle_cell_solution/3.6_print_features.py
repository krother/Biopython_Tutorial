
from Bio.Seq import Seq
from Bio import SeqIO

# task 3.6
records = SeqIO.parse('L26462.gb', "genbank")
globin = records.next()
for f in globin.features:
    print f
