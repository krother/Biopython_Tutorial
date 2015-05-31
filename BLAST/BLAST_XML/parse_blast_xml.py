
# QUESTION:
# what data do the objects (records, hsps etc.) contain?

from Bio.Blast import NCBIXML

xml_file = open("blast_output.xml")
blast_out = NCBIXML.parse(xml_file)

for record in blast_out:
    for alignment in record.alignments:
        for hsp in alignment.hsps:
            # filter by e-value
            if hsp.expect < 0.0001:
                print hsp.score
                print hsp.query
                print hsp.match
                print hsp.sbjct
                print '-' * 70
                print
                
