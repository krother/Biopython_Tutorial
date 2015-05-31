
# Python BLAST Tutorial

## Overview

In this tutorial, you learn how to run BLAST **locally, many times** and to **read BLAST results with Python**. While doing so, you will familiarize with the concept of **program pipelines** and how you can implement them in Python.

The tutorial consists of six parts:

1. Preparations
2. Running local BLAST manually
3. Running local BLAST from Python
4. Running BLAST many times with Python
5. Reading BLAST output with Biopython
6. Plotting the results

----

## Case Study: Plasmodium falciparum

We have the hypothesis that *Plasmodium falciparum* has adapted to the human organism during its long history as a parasite. Specifically, we want to examine whether proteins from *Plasmodium* are more similar to human proteins than one would expect. If this is true we could interpret e.g. that more similar proteins help *Plasmodium* to evade the human immune system.

As a small sample study, we will BLAST a set of peptides from a few *Homo sapiens* proteins against the proteome of *Plasmodium falciparum*. As a control, we will use the proteome of *Schizosaccharomyces pombe*.

----

## 1. Preparations

### 1.1 Check whether BLAST+ is properly installed

Enter the two following commands in a Linux console:

    makeblastdb
    blastp

Both should result in an error message other than *command not found*.

### 1.2 Create a BLAST database for Plasmodium falciparum

Create a BLAST database for the *Plasmodium* proteins. First, open a console and go to the folder *data/* . Type:

    makeblastdb -in Plasmodium_falciparum.fasta -dbtype prot

You should see a message similar to:

    Adding sequences from FASTA; added 5414 sequences in 0.56993 seconds.

### 1.3 Create a BLAST database for the control organism

Please create a BLAST database for *Schizosaccharomyces pombe* as well. 

### 1.4 Questions

* What files have appeared in the *data/* directory?
* Why do we need to create a database first? Why can't BLAST do that right before each query?

----

## 2. Running local BLAST manually

Before running a large series of BLAST experiments, we will run a small sequence as a technical proof of concept. We are using a sequence copied from the *Plasmodium* sequences, so we know that BLAST should generate a 100% match.

### 2.1 Create a query file

Create an empty file *query.seq* in a text editor. Write the following peptide sequence into the file:

    DAAITAALNANAVK

Make sure that there are no other characters in the file (no empty lines or FASTA deflines). Save the file to the *data/* directory.

### 2.2 Running local BLAST against Plasmodium

Go to a console in the *data/* directory and type:

    blastp -query query.seq -db Plasmodium_falciparum.fasta -out output.txt -outfmt 7

### 2.3 Running local BLAST against the control group

Repeat the above query for *Schizosaccharomyces pombe*.

### 2.4 Adjust output formats

Insert different numbers (1-7) for the **outfmt** parameter and re-run the query. 

### 2.5 Questions

* Take a look at the BLAST output. Is the result what you would expect?
* Does the control group support your assumptions so far?
* Which of the output formats do you find the easiest to read?
* Which of the output formats is probably the easiest to read for a program?

----

## 3. Running local BLAST from Python

Now we are going to do exactly the same operation from a Python program. For this we will need the **os module**.

### 3.1 Introduction to the os module

Open the document *pipelines/os_module_puzzle.pdf*. Do the exercise.

### 3.2 Running BLAST from Python

Now we will use the function **os.system** to run BLAST. Create a Python script **run_blast.py** in the *data/* directory. Write the following commands into it:

    import os

    cmd = "blastp -query query.seq -db Plasmodium_falciparum.fasta -out output.txt -outfmt 7"
    os.system(cmd)

Execute the program.

### 3.3 Customizing the query

In order to make the BLAST command in Python more flexible, we will combine it from variables. Change the code to the following:

    db = "Plasmodium_falciparum.fasta"
    cmd = "blastp -query query.seq -db " + db + " -out output.txt -outfmt 7"

### 3.4 More variables

Now add separate variables for the query and the output file name as well.

### 3.5 Additional examples for using os

In the *pipelines/* directory you find more examples using the *os module*. If you like, try them out as well.

### 3.6 Questions

* Is the output of the Python BLAST run identical to the one you did manually? How can you check that?

----

## 4. Running BLAST many times with Python

### 4.1 Creating query files

The file *data/human_peptide.fasta* contains about 2000 peptides. We want to run BLAST for each of them. To do so, we need to write each peptide to a separate file.

First, create a new folder for the query files:

    mkdir data/queries

The Python script **multiblast/split_fasta.py** does that using **Bio.SeqIO**. You can use it by typing in the *multiblast/* directory:

    python split_fasta.py ../data/human_peptide.fasta ../data/queries

If you want, you can try writing that script by yourself.

### 4.2 Validate the queries

Make sure that the query files have been generated and that they are not empty. You can check both with:

    ls -l data/queries
    more data/queries/9568103_99.fasta

### 4.3 Create output directories

Prepare a place where the results from each BLAST run will be stored:

    mkdir data/Plasmodium_out
    mkdir data/Pombe_out

### 4.4 Run BLAST

You can run BLAST for all queries with the script **multiblast/run_blast.py**. It uses **os** for three different things:

1. Reading directory names as command-line parameters
2. Looping through all files in a directory
3. Running the BLAST command

**However, the program is incomplete.**

You need to complete the BLAST command inserting the file names from the given variables. Use the parameter *-outfmt 5* in order to create XML output. We will need this later to read it from Biopython.

When everything is done, you should be able to execute the script with:

    python run_blast.py ../data/queries/ ../data/Plasmodium_falciparum.fasta ../data/Plasmodium_out/

Inspect the result.

----

## 5. Reading BLAST output with Biopython

### 5.1 Reading XML data

Run the program **BLAST_XML/parse_blast_xml.py**.

    python parse_blast_xml.py

### 5.2 Read one of your BLAST result files

Adjust the program to read one of your BLAST output files. Try to figure out how many HSPs there are, and how many are below an e-value threshold of 0.001.

### 5.3 Read all of your BLAST result files

Customize the program to read **all** of your result files. How many hits do you have in total. What is the hit with the highest score?

----

## References

BLAST+ is a new, faster (C++ based) version that replaces BLAST2, as of Oct 2013. Also see: http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download 
