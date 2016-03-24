
from Bio import SeqIO

# task 3.4
num = 0
records = SeqIO.parse("globins.gb", "genbank")
for r in records:
    if 'L26462' in r.annotations['accessions']:
        SeqIO.write(r, 'L26462_{}.gb'.format(num), \
            "genbank")
        num += 1
        print 'FOUND'


