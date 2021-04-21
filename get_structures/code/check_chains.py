#!/usr/bin/env python3
import os
import sys

pdbid=os.environ["id"]
def summarize_file(fhandle, option):

    models = set()
    chains, resids, atoms, hetatm = set(), set(), set(), set()
    has_altloc, has_icode = False, False  # flags only

    model_id = "X   "
    for line in fhandle:
        if line.startswith('ATOM'):
            chains.add(model_id + line[21])

