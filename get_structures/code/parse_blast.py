#!/usr/bin/env python3
import os

#Parse BLAST result
#####
from Bio.Blast import NCBIXML
E_VALUE_THRESH = float(os.environ["ethresh"])
print ("E_VALUE_THRESH = ", E_VALUE_THRESH)
for record in NCBIXML.parse(open("blast_files/blast_results.xml")): 
    if record.alignments: 
       print("\n") 
       print("query: %s" % record.query) 
       for align in record.alignments: 
          for hsp in align.hsps: 
             if hsp.expect < E_VALUE_THRESH:
                print("=====MATCH=====") 
                print("matching pdb|chain (first to appear):", align.title[4:10])
                print("length:", align.length)
                print("e value:", hsp.expect)
                print("query:",hsp.query[0:70],"...")
                print("match:",hsp.match[0:70],"...")
                print("sbjct:",hsp.sbjct[0:70],"...")
                print('\n')
