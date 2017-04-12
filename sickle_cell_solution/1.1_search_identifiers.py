
from Bio import Entrez 

Entrez.email = "krother@academis.eu"

# task 1.1
handle = Entrez.esearch(db="nucleotide", term="sickle cell AND human NOT chromosome", retmax =100)
records = Entrez.read(handle)
print(records['Count'])
identifiers = records['IdList']
print(identifiers)
