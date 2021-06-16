#!/usr/bin/env python3
import os
import sys
#sys.argv[1]
pdbid = os.environ["id"]
pdbf = pdbid+".pdb"
pdbfh = open(pdbf, 'r')
prev_chain = "oogabooga"
chainlist=list()
for line in pdbfh:
    if line.startswith('ATOM'):
        if line[21] != prev_chain:
            chainlist.append(line[21])
            prev_chain=line[21]
pdbfh.close()
print(chainlist)
