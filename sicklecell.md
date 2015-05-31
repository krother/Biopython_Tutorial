
# Project: Diagnosing Sickle Cell Anemia

**Programming tasks with Biopython**

## Goal

Your goal is to develop an experimental test that reveals whether a patient suffers from the hereditary disease sickle cell anemia. The test for diagnosis should use a restriction enzyme on a patients’ DNA sample. For the test to work, you need to know exactly what genetic difference to test against. In this tutorial, you will use Biopython to find out.
The idea is to compare DNA and protein sequences of sickle cell and healthy globin, and to try out different restriction enzymes on them. 

This tutorial consists of four parts:

1. Use the module Bio.Entrez to retrieve DNA and protein sequences from NCBI databases.
2. Use the module Bio.SeqIO to read, write, and filter information in sequence files.
3. Use the modules Bio.Seq and Bio.SeqRecord to extract exons, transcribe and translate them to protein sequences.
4. Use the module re to identify restriction sites.
Have fun!


### What is sickle cell anemia?

At the beginning of the course, watch the 5-minute movie **"Sickle Cell Anemia"** by Paulo César Naoum and Alia F. M. Naoum: **[http://www.youtube.com/watch?v=R4-c3hUhhyc](http://www.youtube.com/watch?v=R4-c3hUhhyc)**.

## 1. Bio.Entrez
**Retrieving DNA and protein sequences**

### 1.1 Search identifiers on NCBI

Search for the cDNA sequence of the sickle cell globin protein from NCBI.
Use the Entrez.esearch function. As keywords, use ‘sickle cell AND human NOT chromosome’. 
Print the resulting database identifiers (not the full sequences).

### 1.2 Retrieve sequences using identifiers

Using the identifiers from task 1.1, retrieve the full sequence with the Entrez.efetch function from the NCBI server. 

Print the identifier and defline for each entry using a for loop. 
The parameter rettype should be ‘fasta’.

### 1.3 Retrieve a single GenBank entry
In the output of task 1.2, locate the cDNA of the sickle cell globin. Use the identifier to download a GenBank entry only for that sequence. Print the entry. The parameter rettype should be ‘gb’.

### 1.4 Write an output file

Save the GenBank entry from task 1.3 to a file ‘sickle.gb’.

### 1.5 Retrieve and write multiple GenBank entries

Retrieve entries for the gene sequences of the human globin family. Find appropriate keywords to limit the search to beta-globin and only complete coding sequences.  to a file.

### Optional tasks for fast programmers: 

- 1.6 Save the retrieved entries to a single FASTA file.
- 1.7 Save each of the beta-globin sequences to a separate GenBank file.
- 1.8 Use Entrez to search and retrieve recent PubMed entries related to malaria and sickle cell anemia.


## 2. Bio.SeqIO

**Reading, writing, and filtering sequence files**

### 2.1 Read a GenBank file

Read a the ‘sickle.gb’ file from task 1.4 using the SeqIO.parse(). The first parameter of parse() is the filename, the second is ‘genbank’. Use the type() and dir() function to find out what the resulting object is and what attributes it has.

### 2.2 Print information for one sequence

Print the id, name and description of the sickle cell globin entry.

### 2.3 Write a FASTA file

Save the GenBank entry from task 2.1 to a FASTA file using the SeqIO.write() function. The first parameter of write() is a list of sequence records, the second a file open for writing, and the third should be ‘fasta’.

### 2.4 Print information for multiple sequences

Print the id, name, and description of all human beta-globins.

### 2.5 Filtering sequence entries

Print the same information as in task 2.4, but do not show non-globin entries : if the description contains either ‘vector’ or ‘isolate’, don’t print anything.

### Optional tasks for fast programmers:

- 2.6 Filter the list of sequence entries even further using your own criteria and save the filtered list to a FASTA file.
- 2.7 Sort the lines in the parse_FASTA_easy/ folder
- 2.8 do the exercise in the parse_FASTA_difficult/ folder


## 3. Bio.Seq and Bio.SeqRecord

**Working with sequences**

### 3.1 The DNA sequence

Print the DNA sequence of the sickle cell globin cDNA. Use the dir() function to find out the name of the attribute.

### 3.2 Transcribe DNA to RNA

Transcribe the sickle cell cDNA sequence to RNA and print it. Use the transcribe() method of a Seq object.

### 3.3 Translate RNA to protein

Translate the sickle cell RNA to a protein sequence and print it. Use the translate() method of a Seq object. Save the protein sequence to a separate file.

### 3.4 Analyze annotations of beta-globin

Find the  beta-globin entry with accession L26462. Use the field r.annotations['accessions'] on a SeqRecord object. Write the entry to a separate GenBank file.

#### Bonus Question:

Print the complete annotations of the record. What data type in Python is it? Identify three methods to use that data type. 

### 3.6 Extract sequence features

Print all features of the L26452 entry. Use the field r.features on a SeqRecord object.

### 3.7 Extract exons

Print all exon features, their attributes start, end and nofuzzy_start, nofuzzy_end. Use the latter as indices to extract portions of the complete sequence. Concatenate all exon sequences to a single string.
#
### 3.8 Display the beta-globin protein sequence

Transcribe and translate the concatenated exon sequence and print it.

### 3.9 Results

Print the sickle cell and healthy beta-globin sequence in subsequent lines or a text file. Also print the two corresponding protein sequences in subsequent lines. What differences do you see?

### Optional tasks for fast programmers:

- 3.10 What are the ‘N’ characters in the sickle cell globin cDNA? How did they affect the transcription/translation?
- 3.11 Do the exercise in the read_alignment/ folder
- 3.12 Figure out how to read an alignment using Biopython
- 3.13 Design a PCR primer that specifically recognizes the sickle cell gene.


## 4. Pattern Matching
**Identification of restriction sites**

### 4.1 Familiarize with regular expression

Do the first 3-4 exercises on the RegexOne website (http://www.regexone.com)

### 4.2 Search for a start codon

Use the re.search() function to locate the start codon (ATG) in the cDNA sequence of healthy beta-globin. The first parameter is a search pattern string, and the second is the string to be searched.

### 4.3 Search for a restriction site

Create a regular expression using the re.compile() function for the restriction enzyme DdeI (cuts at NN^CTNAG). Search with a regular expression in both sickle cell and beta-globin DNA sequences. If the search method returns a match, the match object has a start and a stop attribute. Print the start and stop found in both DNA sequences.

### 4.4 Apply more restriction enzymes

Test patterns for the restriction sites of:

    HinfI (GTNNAC)
    BceAI (ACGGCNNNNNNNNNNNNN)
    BseRI (GAGGAGNNNNNNNNNN)
    EcoRI (GAATTC)
    MstII (CCTNAGG) 

on both DNA sequences. Which restriction enzyme could you use to specifically identify carriers of the sickle cell anemia gene.

### Optional tasks for fast programmers: 

- 4.5 To facilitate the restriction analysis, replace the N's in the sickle cell DNA by the corresponding positions from the healthy DNA. Print the resulting DNA sequence.
- 4.6 Take a look at the other exercise websites given in the regex_links.pdf document.
