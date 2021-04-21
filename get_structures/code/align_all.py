#!/usr/bin/env python3
import os
from Bio import pairwise2, SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord

align = MultipleSeqAlignment([])

#uncomment to use input seq in alignment
#reffile = SeqIO.read(os.environ["reffasta"], format="fasta")
#refseq = Seq(str(reffile.seq))
#align.add_sequence("Input sequence",str(refseq))

with open("pdb_list_align_input.dat", "r") as lines:
    list = lines.read().splitlines()
    for id in list:
        fastaname=id+".fasta"
        newfile = SeqIO.read(fastaname, format="fasta")
        newseq = Seq(str(newfile.seq))
        align.add_sequence(id,str(newseq))

f = open("align_input.fasta", "w")
#f.write(align.format("fasta"))
f.write(format(align, "fasta"))
f.close()

from Bio.Align.Applications import MafftCommandline
mafft_cline = MafftCommandline("mafft", input="align_input.fasta")
stdout, stderr = mafft_cline()
with open("align_output.fasta", "w") as handle:
    handle.write(stdout)
#align = AlignIO.read("align_output.fasta", "fasta")
