
from Bio import Entrez 

Entrez.email = "krother@academis.eu"

# task 1.3
handle = Entrez.efetch(db="nucleotide", \
    id="179408", rettype="gb", \
    retmode="text")
text = handle.read()

# task 1.4
open('sickle.gb', 'w').write(text)
