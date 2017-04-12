from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# task 3.7
records = SeqIO.parse('L26462_exon.fasta', "fasta")
globin = next(records)

rna = globin.seq.transcribe()
print(rna)

protein = rna.translate()
print(protein)
