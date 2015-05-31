import sys
from Bio import SeqIO

def split_fasta (fastafile, dirname):

    infile = open(fastafile)
    for record in SeqIO.parse(infile, "fasta"):
        seq = str(record.seq)
        filename = record.id.strip().replace('|', '_')
        outfile = open ("%s/%s.fasta" % (dirname, filename[2:]), "w") 
        outfile.write(seq)
        outfile.close()    
        

split_fasta(sys.argv[1], sys.argv[2])
