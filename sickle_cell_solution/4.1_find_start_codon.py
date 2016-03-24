
from Bio import SeqIO
import re

# task 4.1
records = SeqIO.parse('L26462.gb', 'genbank')
globin = records.next()

match = re.search('ATG', str(globin.seq))
print match
print match.start(), match.end()

