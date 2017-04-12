
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# task 3.6
records = SeqIO.parse('L26462.gb', "genbank")
entry = next(records)

exon_seq = Seq('')
for f in entry.features:
    if f.type == 'exon':
        print(f.location)
        start = f.location.nofuzzy_start
        end = f.location.nofuzzy_end
        exon_seq += entry.seq[start:end]
print(exon_seq)

rec = [SeqRecord(exon_seq, 'L26462 exon sequence')]
SeqIO.write(rec, 'L26462_exon.fasta', 'fasta')
