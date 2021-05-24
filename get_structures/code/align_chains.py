#!/usr/bin/env python3
import sys
import os
from Bio import pairwise2
from Bio.Seq import Seq
from Bio import SeqIO
import Bio.PDB
pdb_parser = Bio.PDB.PDBParser(QUIET = True)
#import pymolPy3
#pm = pymolPy3.pymolPy3(0)

#load pdbs in 
ref_id=os.environ["pdb_ref_id"]
pdb_l=list()
clean_f = open("pdb_list_clean.dat")
clean = clean_f.read().splitlines()
for line in clean:
    pdb_l.append(line)
clean_f.close()
#print(pdb_l)
ref_pdb=ref_id+".pdb"


#align structures to ref
#get seq for ref 
#make res sets for ref
ref_res_set=set()
ref_chn_d = {}
ref_chn_l=list()
for chn_id in sys.argv[1:]:
    ref_chn_l.append(chn_id)
##get atoms to align in ref pdb 
#ref_f = open(ref_pdb, 'r')
#for line in ref_f:
#    if line.startswith('ATOM'):
#        #all res
#        ref_res_set.add(line[22:26])
#        #each chain
#        for chn_id in ref_chn_l:
#            if line[21].strip()== chn_id:
#                ref_chn_d[chn_id].add(line[22:26])
#                if line[12:16].strip()== 'CA':
#                    atoms_chm_d[chn_id].add()
#ref_f.close()

#make dictionary containing the lists of all CA atoms for each chain
ref_structure = pdb_parser.get_structure("reference", ref_pdb)
ref_model    = ref_structure[0]
ref_res_d = {}
for chn_id in ref_chn_l:
    ref_res_d["{0}".format(chn_id)]=[]
    for ref_chain in ref_model:
        if ref_chain.get_id().strip() == chn_id:
            for ref_res in ref_chain:
                try:
                    ref_res_d["{0}".format(chn_id)].append(ref_res['CA'])
                except:
                    continue

#get seqs for ref
ref_records=SeqIO.parse(ref_pdb, "pdb-atom")

#align chains
for pdb in pdb_l:
    pdbf=pdb+".pdb"
    sample_records=SeqIO.parse(pdbf, "pdb-atom")
    for chn_id in ref_chn_l:
        chain_to_be_aligned=chn_id
        sample_structure = pdb_parser.get_structure("sample", pdbf)
        
        # Use the first model in the pdb-files for alignment for formatting (all pdbs should have 1 model anyway)
        sample_model = sample_structure[0]
           
        sample_atoms = []
        for sample_chain in sample_model:
            if sample_chain.get_id().strip() == chn_id:
                for sample_res in sample_chain:
                    try:
                        sample_atoms.append(ref_res['CA'])
                    except:
                        continue

        #pairwise seq align to trim based on 
        for record in ref_records:
            if record.annotations["chain"] == chn_id:
                ref_seq = record.seq
        for record in sample_records:
            if record.annotations["chain"] == chn_id:
                sample_seq = record.seq

            #local finds best matching subseq
            alignments = pairwise2.align.localxx(ref_seq, sample_seq)
            for alignment in alignments:
                print(pairwise2.format_alignment(*alignments[0]))

        #trim residues of sample list to match number of atoms in ref list
        #or should i align structures based on paired residues? 

        #align structures 
        super_imposer = Bio.PDB.Superimposer()
        super_imposer.set_atoms(ref_atoms, sample_atoms)
        super_imposer.apply(sample_model.get_atoms())
        
        # Print RMSD:
        print(super_imposer.rms)
        
        # Save the aligned version
        io = Bio.PDB.PDBIO()
        io.set_structure(sample_structure) 
        io.save(pdb+"_aligned.pdb")	

quit()
#get seq for each chain

#compare each ref chain to each chain of other pdbs
	#seq 
	#distance com
	#number of residues
	#align struc
#rename chains to be the same name as the ref chain it was most similar 

#

