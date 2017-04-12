
from Bio.Seq import Seq
import re

# task 4.3

globin = 'ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTG'
sickle = 'ATGGTGCATCTGACTCCTGTGGAGAAGTCYGCYGTTACTG'
#sickle = 'ATGGTNCAYYTNACNCCNGTGGAGAAGTCYGCYGTNACNG'

HinfI = re.compile('GT..AC')
BceAI = re.compile('ACGGC..............')
BseRI = re.compile('GAGGAG..........')
EcoRI = re.compile('GAATTC')
MstII = re.compile('CCT.AGG')

ENZYMES = [HinfI, BceAI, BseRI, EcoRI, MstII]

for enz in ENZYMES:
    r1 = enz.search(globin)
    r2 = enz.search(sickle)
    if r1:
        print('GLOBIN:', r1.start(), r1.end())
    if r2:
        print('SICKLE:', r2.start(), r2.end())
