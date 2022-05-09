
# Biopython Examples

## 1. Getting started

    :::python
    import Bio
    from Bio.Seq import Seq
    dna = Seq("ACGTTGCAC")
    print(dna)

#### (alternative)

    :::python
    from Bio.Alphabet import IUPAC
    dna = Seq("AGTACACTGGT", IUPAC.unambiguous_dna)


## 2. Reverse complement, transcribing & translating

    :::python
    dna.reverse_complement()
    rna = dna.transcribe()
    rna.translate()


#### (alternative)

    :::python
    from Bio.Seq import reverse_complement, transcribe, translate
    reverse_complement("GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCT")

## 3. Calculating GC-content

    :::python
    from Bio.SeqUtils import GC
    GC(dna)

## 4. Caculating molecular weight (DNA only)

    :::python
    from Bio.SeqUtils import molecular_weight
    molecular_weight("ACCCGT")

## 5. Loading sequences from a FASTA file

    :::python
    from Bio import SeqIO
    for record in SeqIO.parse("ls_orchid.fasta", "fasta"):
        print record.seq, len(record.seq)

## 6. Plotting a histogram of seq lengths with pylab 

`pylab` aka `matplotlib` needs to be installed separately.

    :::python
    import pylab
    sizes=[len(r.seq) for r in SeqIO.parse("ls_orchid.fasta","fasta")]
    pylab.hist(sizes, bins=20)
    pylab.title("%i orchid sequences\nLengths %i to %i" \
                % (len(sizes), min(sizes), max(sizes)))
    pylab.xlabel("Sequence length (bp)")
    pylab.ylabel("Count")
    pylab.show()
