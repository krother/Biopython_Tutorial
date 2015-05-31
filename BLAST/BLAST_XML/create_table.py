
import os
import sys

# uses data created with -outfmt 6

def create_table (dir1, dir2):
    
    filelist1 = os.listdir(dir1)
    filelist1.sort()
    filelist2 = os.listdir(dir2)
    filelist2.sort()
    peptides = []
    plasmodium = []
    pombe = []
    
    for fname in filelist1:
        pepname = fname.split('.')[0]
        peptides.append(pepname)

    for fname1 in filelist1:
        firstline = open(dir1+fname1).readline()
        if firstline.strip() == '':
            score = 0.0
            plasmodium.append(score)
        else:
            columns = firstline.split()
            score = float (columns[-1])
            plasmodium.append(score)

    for fname2 in filelist2:
        firstline = open(dir2+fname2).readline()
        if firstline.strip() == '':
            score = 0.0
            pombe.append(score)
        else:
            columns = firstline.split()
            score = float (columns[-1])
            pombe.append(score)

    print peptides, plasmodium, pombe

create_table("Plasmodium_out/","Pombe_out/")
