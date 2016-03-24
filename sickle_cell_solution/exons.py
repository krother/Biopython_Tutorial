
# extract exons from L26462.gb

import Bio
from Bio import SeqIO

gen = SeqIO.parse('L26462.gb', 'genbank')

entry = gen.next()
sequence = entry.seq

result = ''

for feature in entry.features:
    if 'exon' in feature.type:
        print feature.location
        start = feature.location.nofuzzy_start
        end = feature.location.nofuzzy_end
        result = result + sequence[start:end]
    
print result
print result.transcribe().translate()

from Bio.Seq import Seq

sickle = Seq('ATGGTNCAYYTNACNCCNGTGGAGAAGTCYGCYGTNACNGCNCTNTGGGGYAAGGTNAAYGTGGATGAAGYYGGYGGYGAGGCCCTGGGCAGNCTGCTNGTGGTCTACCCTTGGACCCAGAGGTTCTTNGANTCNTTYGGGGATCTGNNNACNCCNGANGCAGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTNCAYGTGGATCCTGAGAACTTCAGGCTNCTNGGCAACGTGYTNGTCTGYGTGCTGGCCCATCACTTTGGCAAAGAATTCACCCCACCAGTGCANGCNGCCTATCAGAAAGTGGTNGCTGGTGTNGCTAATGCCCTGGCCCACAAGTATCACTAAGCTNGCYTTYTTGYTGTCCAATTT')
print sickle.translate()






