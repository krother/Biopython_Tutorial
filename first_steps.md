
# First Steps in Biopython

## Preparations

1. Start a Python console
2. Type `import Bio`

----

## 1. Load a sequence file

Load the FASTA file :::file data/ap006852.fasta into Biopython.

    :::python
    from Bio import SeqIO

    records = list(SeqIO.parse("ap006852.fasta", "fasta"))
    dna = records[0]

    print(dna.name)
    print(dna.description)
    print(dna.seq[:100])

Check whether the following statements are `True` or `False`:

* The command `print(len(dna))` displays the length of the sequence.
* Replacing `records[0]` by `records[1]` results in a different sequence record.
* Replacing `dna.seq[:100]` by `dna.seq[50:100]` displays a different portion of the sequence.

----

## 2. Manipulating a sequence

    :::python
    from Bio.Seq import Seq
    dna = Seq("ACGTTGCAC")
    print(dna)

    result = dna.__________()
    print(result)
    if result == 'GTGCAACGT':
        print('OK')

**Which of the following needs to be inserted to obtain `GTGCAACGT`?**

1. reverse_complement
2. transcribe
3. translate

----

## 3. Calculating GC-content

Look up *Section 3.2* of the Biopython documentation on ([http://biopython.org/DIST/docs/tutorial/Tutorial.html](http://biopython.org/DIST/docs/tutorial/Tutorial.html)) to find out how to calculate the GC-content of a sequence.

What is the GC-content of the sequence loaded in task 1?

* 55.556
* 46.875
* 33.514
* 50.000

----

## 4. Print annotation of a GenBank file

Load the GenBank file :::file data/ap006852.gbk . In contrast to a FastA file, this one contains not only the sequence, but a rich set of annotations. Load the file as follows:

    :::python
    records = list(SeqIO.parse("ap006852.gbk", "genbank"))
    dna = records[0]

Answer the following questions:

#### 4.1 Which command gives the species the sequence is from?

    :::python
    print(dna.annotations['species'])

    print(dna.annotations['organism'])

#### 4.2 Which command produces the bigger ID number (GenBank ID or PubMed ID)

    :::python
    print(dna.annotations['gi'])

    print(dna.annotations['references'][0])

#### 4.3 How can you view a list of available annotation fields?

    :::python
    print(dna.annotations.get('keys'))

    print(dna.annotations.keys())

----

## 5. Count atoms in a PDB structure

The following code reads the 3D structure of a tRNA molecule from the file :::file data/1ehz.pdb and counts the number of atoms.

**There is a bug in the program. Execute the program. Identify the problem and fix it.**

    :::python
    from Bio import PDB

    parser = PDB.PDBParser()
    struc = parser.get_structure("tRNA", "1ehz.pdb")

    n_atoms = 0
    for model in struc:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    print(residue.resname, residue.id)
                n_atoms += 1

    print(n_atoms)

----

## 6. Retrieve entries from NCBI databases

Use the following code to download identifiers (with the `esearch` web app) and protein sequences for these identifiers (with the `efetch` web app) from the NCBI databases.

**The order of lines got messed up! Please sort the lines to make the code work.**

    :::python
    identifiers = records['IdList']

    from Bio import Entrez

    records = handle.read()

    handle = Entrez.efetch(db="protein", id=identifiers, rettype="fasta", retmode="text")

    open("globins.fasta", "w").write(records)

    searchresult = Entrez.esearch(db="protein", term="hemoglobin", retmax=5)

    records = Entrez.read(searchresult)

    Entrez.email = "krother@academis.eu"

