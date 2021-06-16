#!/usr/bin/env python3
import os
import sys
import requests
import re 
from urllib.error import HTTPError

pdbid=os.environ["id"]

if not re.match(r'[0-9a-zA-Z]{4}$', pdbid):
    emsg = 'invalid pdb code: \'{}\'\n'
    sys.stderr.write(emsg.format(pdbid))
    sys.exit(1)

base_url = 'https://files.rcsb.org/download/'
pdb_type = '.pdb'
pdb_url = base_url + pdbid.lower() + pdb_type

filename=os.environ["id"] + '.pdb'
r=requests.get(pdb_url)
if r.status_code == 404:
    emsg = 'pdb code not found: \'{}\'\n'
    sys.stderr.write(emsg.format(pdbid))
    sys.exit(1)

open(filename, 'wb').write(r.content)

