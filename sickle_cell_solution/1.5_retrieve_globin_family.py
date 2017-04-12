
from Bio import Entrez 

Entrez.email = "krother@academis.eu"

# task 1.5

# search identifiers
handle = Entrez.esearch(db="nucleotide", term="beta-globin AND human AND complete cds NOT chromosome NOT thalassemia NOT zebra fish", retmax=150)
records = Entrez.read(handle)
print(records['Count'])
identifiers = records['IdList']
print(identifiers)

# retrieve sequence entries
handle = Entrez.efetch(db="nucleotide", id=identifiers, rettype="gb", retmode="text")
records = handle.read()
open("globins.gb", "w").write(records)
