#!/usr/bin/env python3
import os

#Send seq to BLAST remotely
#####
from Bio.Blast import NCBIWWW
from Bio import SeqIO 
record = SeqIO.read(os.environ["fasta"], format="fasta")
result_handle = NCBIWWW.qblast("blastp", "pdb", record.format("fasta"), expect=float(os.environ["ethresh"]))
with open("blast_files/blast_results.xml", "w") as out_handle: 
    out_handle.write(result_handle.read())
result_handle.close() 
