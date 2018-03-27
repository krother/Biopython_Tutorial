
# Using NCBI E-utilities

## Using Entrez from Biopython

### Step 1: import Entrez

    from Bio import Entrez 

### Step 2: enter your e-mail

The NCBI server might block anonymous requests, especially big ones!

    Entrez.email = "my@email.eu"


### Step 3: Call esearch to find IDs

    handle = Entrez.esearch(db="value", term="keywords", retmax=100)

Parameters include:

| parameter | examples |
|-----------|----------|
| db        | nucleotide |
|           | protein    |
|           | pubmed     |
| term      | human[Organism] |
|           | hemoglobin |
|           | hemoglobin AND alpha |
| retmax    | 10 (identifiers returned) |

### Step 4: get a list of IDs out of esearch

    records = Entrez.read(handle)
    identifiers = records['IdList']

### Step 5: use efetch to retrieve entries

We use the list of identifiers from step 4:

    handle = Entrez.efetch(db="value", id=identifiers, retmax="200", 
             rettype="fasta", retmode="text")

To read data from text entries as a string:

    text = handle.read()

To read records from XML entries:

    records = Entrez.read(handle)

In addition to the above, parameters include:

| parameter | examples |
|-----------|----------|
| id        | single id |
| rettype   | fasta    |
|           | gb       |
| retmode   | text     |
|           | xml      |


## Documentation:

You find a full list of available options on 
[http://www.ncbi.nlm.nih.gov/books/NBK25500/](http://www.ncbi.nlm.nih.gov/books/NBK25500/)

## Example URLs

### 1. Searching for papers in PubMed

[http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=thermophilic,packing&rettype=uilist](http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=thermophilic,packing&rettype=uilist)
    
### 2. Retrieving publication records in Medline format

[http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=11748933,11700088&retmode=text&rettype=medline](http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=11748933,11700088&retmode=text&rettype=medline)

### 3. Searching for protein database entries by keywords

[http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term=cancer+AND+human](http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term=cancer+AND+human)

### 4. Retrieving protein database entries in FASTA format

[http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=1234567&rettype=fasta](http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=1234567&rettype=fasta)

### 5. Retrieving protein database entries in Genbank format

[http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=1234567&rettype=gb](http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=1234567&rettype=gb)

### 6. Retrieving nucleotide database entries

[http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=9790228&rettype=gb](http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=9790228&rettype=gb)

