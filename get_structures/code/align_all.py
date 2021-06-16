#!/usr/bin/env python3
import os
import sys
from Bio import pairwise2, SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord

align = MultipleSeqAlignment([])


#uncomment to use input seq in alignment
#reffile = SeqIO.read(os.environ["reffasta"], format="fasta")
#refseq = Seq(str(reffile.seq))
#align.add_sequence("Input sequence",str(refseq))

pdblistfile=sys.argv[1]
with open(pdblistfile, "r") as lines:
    list = lines.read().splitlines()
    for id in list:
        fastaname=id+".fasta"
        newfile = SeqIO.read(fastaname, format="fasta")
        newseq = Seq(str(newfile.seq))
        align.add_sequence(id,str(newseq))

f = open("align_input.fasta", "w")
f.write(format(align, "fasta"))
f.close()

#align settings
from Bio.Align.Applications import MafftCommandline
mafft_cline = MafftCommandline("mafft", input="align_input.fasta")
mafft_cline.op = 1.53
mafft_cline.ep = 0.123
mafft_cline.maxiterate = 10 

#do alignment
stdout, stderr = mafft_cline()
with open("align_output.fasta", "w") as handle:
    handle.write(stdout)
#align = AlignIO.read("align_output.fasta", "fasta")
