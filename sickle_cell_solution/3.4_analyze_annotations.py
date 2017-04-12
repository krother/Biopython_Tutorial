
from Bio import SeqIO

# task 3.4
num = 0
records = SeqIO.parse("globins.gb", "genbank")
for rec in records:
    if 'L26462' in rec.annotations['accessions']:
        print(rec)
        SeqIO.write(rec, 'L26462_{}.gb'.format(num), "genbank")
        num += 1
