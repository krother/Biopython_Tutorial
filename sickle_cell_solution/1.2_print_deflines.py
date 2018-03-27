
from Bio import Entrez 

Entrez.email = "krother@academis.eu"

# task 1.1
handle = Entrez.esearch(db="nucleotide", term="human[Organism] AND sickle cell AND globin", retmax =100)
records = Entrez.read(handle)
identifiers = records['IdList']

# task 1.2
handle = Entrez.efetch(db="nucleotide", id=identifiers, retmax="200", rettype="fasta", retmode="xml")
records = Entrez.read(handle)
rec = list(records)
print(rec[0].keys())
for r in rec:
    print(r['TSeq_accver'], r['TSeq_defline'])
