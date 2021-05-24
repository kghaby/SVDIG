#!/usr/bin/env python3
import os
from Bio.PDB import *
from Bio.Blast import NCBIXML

pdb_list=set()
E_VALUE_THRESH = float(os.environ["newethresh"])
for record in NCBIXML.parse(open("blast_files/blast_results.xml")):
    if record.alignments:
       for align in record.alignments:
          for hsp in align.hsps:
             if hsp.expect < E_VALUE_THRESH:
                pdb_list.add(align.title[4:8])

for x in pdb_list:
    print(x)
