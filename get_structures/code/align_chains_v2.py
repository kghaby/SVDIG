#!/usr/bin/env python3
import os
#os.chdir("")
import sys
import pandas as pd
from Bio import pairwise2
from Bio.Seq import Seq
from Bio import SeqIO
import Bio.PDB
import warnings
from Bio.PDB.PDBParser import PDBConstructionWarning
from Bio.SeqIO.PdbIO import BiopythonParserWarning
warnings.simplefilter('ignore', PDBConstructionWarning)
warnings.simplefilter('ignore', BiopythonParserWarning)
pdb_parser = Bio.PDB.PDBParser(QUIET = True)
super_imposer = Bio.PDB.Superimposer()
from Bio import Align
aligner = Align.PairwiseAligner()
from Bio.Align import substitution_matrices

#align params based on Smith-Waterman local alignment alg
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5
aligner.target_end_gap_score = 0.0
aligner.query_end_gap_score = 0.0
aligner.mode = 'local'
#aligner.match_score = 1 #auto
#aligner.mismatch_score = 0 #auto

def find_res(res_nbr,given_list):
    found=False
    for y in given_list:
        for x in y:
            for i in x:
                pos=x.index(i)
                if (pos % 2) == 0:
                    i1=i
                else:
                    i2=i
                    if found == False:
                        found=res_nbr in range(i1,i2)
                        if found:
                            return(found)

def min_rmsd(i):
    ref_atoms=[]
    for ref_chain in ref_model:
        if ref_chain.get_id().strip() == df.ref[i]:
            for ref_res in ref_chain:
                for record in ref_records: #get residue start
                    if record.annotations["chain"] == df.ref[i]:
                        r_res=(ref_res.id[1]-record.annotations["start"])
                r_res_matched_l=r_res_matched_d["{0}".format(str(df.ref[i]+"_"+df.samp[i]))]
                if find_res(r_res,r_res_matched_l): #if the res is one that was aligned (not a gap)
                    try:
                        ref_atoms.append(ref_res['CA'])
                    except:
                        continue
    sample_atoms = []
    for sample_chain in sample_model:
        if sample_chain.get_id().strip() == df.samp[i]:
            for sample_res in sample_chain:
                for record in sample_records: #get residue start
                    if record.annotations["chain"] == df.samp[i]:
                        s_res=(sample_res.id[1]-record.annotations["start"])
                s_res_matched_l=s_res_matched_d["{0}".format(str(df.ref[i]+"_"+df.samp[i]))]
                if find_res(s_res,s_res_matched_l):
                    try:
                        sample_atoms.append(sample_res['CA'])
                    except:
                        continue

    super_imposer.set_atoms(ref_atoms, sample_atoms)
    super_imposer.apply(sample_model.get_atoms())
    return(super_imposer.rms)

def rename_chains(pdbfh, chains):
    chain_from, chain_to = chains
    #records = ('ATOM', 'HETATM', 'TER')
    ter = True
    for line in pdbfh:
        if line.startswith('ATOM'):
            if line[21] == chain_from:
                yield line[:21] + chain_to + line[22:]
                ter = False
                continue
        elif line.startswith('TER'):
            if not ter:
                yield 'TER\n'
                ter = True
                continue
            else:
                continue
        else:
            continue


#load pdbs in
ref_id=os.environ["pdb_ref_id"] #ref_id="3JA8" #os.thing
pdb_l=list()
clean_f = open("pdb_list_standard_input.dat")
clean = clean_f.read().splitlines()
for line in clean:
    pdb_l.append(line)
clean_f.close()

ref_pdb=ref_id+".pdb"
ref_structure = pdb_parser.get_structure("reference", ref_pdb)
ref_model = ref_structure[0]

#align chains
for pdb in pdb_l:
    print("=================")
    print(pdb,"-",ref_id)
    nomore=False
    #get ref chains (at start of every loop bc its slowly deleted during loop)
    ref_chn_l=list()
    ref_records=list(SeqIO.parse(ref_pdb, "pdb-atom"))
    for record in ref_records:
        ref_chn_l.append(record.annotations["chain"])

    pdbf=pdb+".pdb"
    sample_records=list(SeqIO.parse(pdbf, "pdb-atom"))
    sample_structure = pdb_parser.get_structure("sample", pdbf)
    # Use the first model in the pdb-files for alignment for formatting (all pdbs should have 1 model anyway)
    sample_model = sample_structure[0]
    s_chn_l=list()
    ref_matched=''
    ref_matched_l=[]
    r_res_matched_d={}
    s_matched=''
    s_matched_l=[]
    s_res_matched_d={}
    df_best=pd.DataFrame(columns=['pdb', 'samp', 'ref', 'score'])
    for s_record in sample_records:
        s_chn_l.append(s_record.annotations["chain"])
    l_r=len(ref_chn_l)
    for nbr_chains in range(l_r):
        df=pd.DataFrame(columns=['pdb', 'samp', 'ref', 'score'])
        if len(s_chn_l) == 0: #if used all chains in samp already, move on to the next pdb
            nomore=True
            break
        for chn_id in ref_chn_l:
            bestscore=0
            #pairwise seq align
            for r_record in ref_records:
                if r_record.annotations["chain"] in ref_matched_l: #if record belongs to a chain that was already matched, skip
                    continue
                if r_record.annotations["chain"] == chn_id:
                    ref_seq = r_record.seq
            for s_record in sample_records:
                if s_record.annotations["chain"] in s_matched_l: #if record belongs to a chain that was already matched, skip
                    continue
                sample_seq = s_record.seq
                #local finds best matching subseq
                #alignments = pairwise2.align.localxx(ref_seq, sample_seq) #bad aligner
                alignments = aligner.align(ref_seq, sample_seq)
                r_res_matched_l=[alignments[0].aligned[0]]
                s_res_matched_l=[alignments[0].aligned[1]]
                r_res_matched_d["{0}".format(str(chn_id+"_"+s_record.annotations["chain"]))]=r_res_matched_l
                s_res_matched_d["{0}".format(str(chn_id+"_"+s_record.annotations["chain"]))]=s_res_matched_l
                score=alignments[0].score
                if score > bestscore:
                    bestscore=score #newbest
                    bestchain=s_record.annotations["chain"].strip()
                if s_record.annotations["chain"] == s_chn_l[-1]: #if last chain, append best match to df
                    df=df.append({'pdb':pdb, 'samp':bestchain, 'ref':chn_id, 'score':bestscore},ignore_index=True) #best sample-ref matches for all (remaining) chains. may have dup samp chains
        #get best match of whole matches to enforce bijection
        matches=0
        for i in range(len(df)):
            if df.score[i] == df.max().score:
                matches=matches+1
                if matches > 1: #if more than one chain matches...
                    #align structures to get rmsd to break tie
                    print("Breaking tie between",pdb,df.samp[i],"-",ref_id,df.ref[i],"and",pdb,df.samp[prev_i],"-",ref_id,df.ref[prev_i],"using RMSDs")
                    i_rmsd=min_rmsd(i)
                    prev_i_rmsd=min_rmsd(prev_i)
                    print("\t",pdb,df.samp[i],"-",ref_id,df.ref[i],":",i_rmsd)
                    print("\t",pdb,df.samp[prev_i],"-",ref_id,df.ref[prev_i],":",prev_i_rmsd)
                    if i_rmsd < prev_i_rmsd: #smaller RMSD => better match
                        print('\t\t',pdb,df.samp[i],"-",ref_id,df.ref[i],'is the better match, so these chains are getting credit first.')
                        s_matched=df.samp[i]
                        ref_matched=df.ref[i]
                        sc=df.score[i]
                        prev_i=i #last i to match
                        continue
                    elif i_rmsd > prev_i_rmsd:
                        print('\t',pdb,df.samp[prev_i],"-",ref_id,df.ref[prev_i],"is the better match")
                        continue #all variables already set from the first matching loop
                    elif i_rmsd == prev_i_rmsd:
                        print("WARNING: cant tie break due to equal RMSDs. using first match")
                        print('\t',"i-row:",df[i])
                        print('\t',"prev_i-row:",df[prev_i])
                        print("i didnt think this would happen...")
                        print("if you want to add another tiebreaker, add another if statement to code/align_chains_v2.py under '#other tie breakers here'")
                        #other tie breakers here
                        continue
                s_matched=df.samp[i]
                ref_matched=df.ref[i]
                sc=df.score[i]
                prev_i=i #last i to match
        df_best=df_best.append({'pdb':pdb, 'samp':s_matched, 'ref':ref_matched, 'score':sc},ignore_index=True)
        #remove matched chains before next iteration
        ref_chn_l.remove(ref_matched)
        s_chn_l.remove(s_matched)
        ref_matched_l.append(ref_matched)
        s_matched_l.append(s_matched)
    print(df_best)
    if nomore:
        print("no more chains in sample pdb",pdb,"to compare to ref pdb",ref_id)
    if len(df_best) > l_r:
        print("WARNING: df_best (ie the printed table ie best of the best matches) for",pdb,"contains more rows(",len(df_best),") than chains in reference (",l_r,"), implying nonbijective matching. Skipping")
        continue

    #rename chains based on df_best and only keep matched chains
    print('\nchanging chain names...')
    pdb_stndrd = open("{}_stndrd.pdb".format(pdb), "a")

    ref_chains_ordered = sorted(ref_matched_l)
    for chain_to in ref_chains_ordered: #each match in alphabetical/numerical order
        pdbfh = open(pdbf, 'r')
        try:
            chain_from=df_best.loc[df_best['ref'] == chain_to, 'samp'].iloc[0]
        except:
            continue
        chains = (chain_from,chain_to)
        print('\t',"writing {0} as {1}".format(chain_from,chain_to))
        renamed_chain = rename_chains(pdbfh, chains)
        for line in renamed_chain:
            pdb_stndrd.write(line)
        pdbfh.close()
    pdb_stndrd.close()
