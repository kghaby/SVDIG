#!/usr/bin/env python3
import os
from Bio import pairwise2
from Bio.Seq import Seq
from Bio import SeqIO

#pairwise to input seq
reffile = SeqIO.read(os.environ["reffasta"], format="fasta")
newfile = SeqIO.read(os.environ["newfasta"], format="fasta")
refseq = Seq(reffile) 
newseq = Seq(newfile)
alignments = pairwise2.align.globalxx(refseq, newseq)
#test_alignments = pairwise2.align.localds(refseq, newseq, blosum62, -10, -1)
for alignment in alignments: 
    print(format_alignment(*alignment))

