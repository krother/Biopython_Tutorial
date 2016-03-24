
from Bio import Entrez 

Entrez.email = "krother@academis.eu"

# task 1.3
handle = Entrez.efetch(db="nucleotide", id="179408", rettype="gb", retmode="text")
text = handle.read()
print text
