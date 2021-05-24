#!/usr/bin/env python3
import os

#Parse BLAST result
#####
from Bio.Blast import NCBIXML
E_VALUE_THRESH = float(os.environ["ethresh"])
for record in NCBIXML.parse(open("blast_files/blast_results.xml")): 
    if record.alignments: 
       print("PDB|CHAIN","\t","Expectation value")
       for align in record.alignments: 
          for hsp in align.hsps: 
             print(align.title[4:10],"\t",hsp.expect)
