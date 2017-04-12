from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# task 3.8    
records = SeqIO.parse('L26462_exon.fasta', "fasta")
globin = records.next()

rna = globin.seq.transcribe()
print rna

protein = rna.translate()
print protein
