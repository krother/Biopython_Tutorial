
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# task 3.9
records = SeqIO.parse('L26462_exon.fasta', "fasta")
globin = records.next()

records = SeqIO.parse('sickle.gb', "genbank")
sickle = records.next()

print globin.seq[:70]
print sickle.seq[:70]

print

prot1 = globin.seq.transcribe().translate()
prot2 = sickle.seq.transcribe().translate()

print prot1[:70]
print prot2[:70]
