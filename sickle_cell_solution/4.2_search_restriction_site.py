
from Bio.Seq import Seq
import re

# task 4.2
DdeI = re.compile('CT.AG')

globin = 'ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTG'
sickle = 'ATGGTGCATCTGACTCCTGTGGAGAAGTCYGCYGTTACTG'

r1 = DdeI.search(globin)
r2 = DdeI.search(sickle)
if r1:
    print('GLOBIN', r1.start(), r1.end())
if r2:
    print('SICKLE', r2.start(), r2.end())
