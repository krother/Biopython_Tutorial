
from Bio.Seq import Seq
from Bio import SeqIO

# task 3.5
records = SeqIO.parse('L26462.gb', "genbank")
globin = next(records)
for f in globin.features:
    print(f)
