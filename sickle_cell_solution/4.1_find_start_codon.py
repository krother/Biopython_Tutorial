
from Bio import SeqIO
import re

# task 4.1
records = SeqIO.parse('L26462.gb', 'genbank')
globin = next(records)

match = re.search('ATG', str(globin.seq))
if match:
    print(match.start(), match.end())
