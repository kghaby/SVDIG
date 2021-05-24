import os
import sys
with open(os.environ["chnfile"], "r") as fhandle:
    records = ('ATOM')
    sequence = []  # list of chain sequences
    seen = set()
    prev_chain = None
    for line in fhandle:
        if line.startswith(records):
            chain_id = line[21]
            if chain_id != prev_chain:
                sequence.append([chain_id])
                prev_chain = chain_id

    labels = sorted(set([c[0] for c in sequence]))           
    print(''.join(labels))
