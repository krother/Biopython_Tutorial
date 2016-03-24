
from Bio import SeqIO

# task 2.5
result = []
records = SeqIO.parse("globins.gb", "genbank")
for r in records:
    if not ('vector' in r.description or 'isolate' in r.description):
        print 'ACCEPT', r.description	
        result.append(r)
    else:
        print 'SKIP', r.description	

print len(result)
out = open('globins_filtered.fasta', 'w')
SeqIO.write(result, out, "fasta")	
