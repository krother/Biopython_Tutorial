
# Biopython – Basics

Biopython is a Python library for reading and writing many common biological data formats. It contains some functionality to perform calculations, in particular on 3D structures. The library can be found at [www.biopython.org](www.biopython.org).

## 1. Getting started

    import Bio
    from Bio.Seq import Seq
    dna = Seq("ACGTTGCA")
    print(dna)

#### (alternative)

    from Bio.Alphabet import IUPAC
    dna = Seq("AGTACACTGGT", IUPAC.unambiguous_dna)


## 2. Reverse complement, transcribing & translating

    dna.reverse_complement()
    rna = dna.transcribe()
    rna.translate()


#### (alternative)

    from Bio.Seq import reverse_complement, transcribe, translate
    reverse_complement("GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCT")

## 3. Calculating GC-content

    from Bio.SeqUtils import GC
    GC(dna)

## 4. Caculating molecular weight (DNA only)

    from Bio.SeqUtils import molecular_weight
    molecular_weight("ACCCGT")

## 5. Loading sequences from a FASTA file

    from Bio import SeqIO
    for record in SeqIO.parse("ls_orchid.fasta", "fasta"):
        print record.seq, len(record.seq)

## 6. Plotting a histogram of seq lengths with pylab (needs to be installed separately)

    import pylab
    sizes=[len(r.seq) for r in SeqIO.parse("ls_orchid.fasta","fasta")]
    pylab.hist(sizes, bins=20)
    pylab.title("%i orchid sequences\nLengths %i to %i" \
                % (len(sizes), min(sizes), max(sizes)))
    pylab.xlabel("Sequence length (bp)")
    pylab.ylabel("Count")
    pylab.show()
