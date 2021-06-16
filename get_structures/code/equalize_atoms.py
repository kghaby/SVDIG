#!/usr/bin/env python3
import os
import sys
import re
from Bio import AlignIO

af = sys.argv[1] #"../mcm_test/pdb_bank/aligned/alignment_files/align_output.fasta" #sys.argv[1]
align = AlignIO.read(af, "fasta") 
pdbid= os.environ["id"] #"3JA8_stndrd" #os.environ["id"]
pdb=pdbid + ".pdb"
prev_length="oogabooga"

atomlist = sys.argv[2:] #["N","CA","C","O"] #sys.argv[2:]

#get alignment length and check to make sure its the same for each row
for record in align:
    if prev_length != "oogabooga":
        if len(record.seq) != prev_length:
            emsg = 'ERROR: rows in alignment output file differ in length\n'
            sys.stderr.write(emsg)
            sys.exit(1)
    if record.id == pdbid:
        alignment_length=len(record.seq)
    prev_length=len(record.seq)
    
if alignment_length > 9999:
            emsg = 'ERROR: alignment length is greater than 9999, so it cannot fit in 1 pdb\n'
            sys.stderr.write(emsg)
            sys.exit(1)

#write new pdb         
newpdb=pdbid + "_withgaps" + ".pdb"
f = open(newpdb, "w")
pdbfh = open(pdb, 'r')
ter=False
terline='TER'+' '*77
gapchain=1
gapres=-1
gapatom='X'
prev_chain=None
col_index=1
prev_res=99999999
#gapline='ATOM      0  X   GAP 0   0       nan     nan     nan    nan   nan'
def gapline_func(gapchain,gapres,gapatom):
    gapline='ATOM      0  '+str(gapatom).ljust(3)+' GAP '+str(gapchain)+str(gapres).rjust(4)+'       nan     nan     nan    nan   nan'
    return gapline
gapped=False
for record in align:
    if record.id == pdbid:
        for line in pdbfh:
            if line.startswith('ATOM'):
                if line[12:16].strip() not in atomlist:
                    continue
                if prev_res==int(line[22:26]): #if same res, write line and avoid iterating through alignment columns
                    f.write(line)
                    ter=False
                    gapped=False
                    continue
                for n in range(col_index, alignment_length+1):
                    if record.seq[n-1:n] == '-':
                        #write gapline
                        if prev_chain == None: #name chain for getchain 
                            gapchain=line[21]
                        elif prev_chain != line[21]:
                            gapchain=prev_chain
                        elif prev_chain == line[21]:
                            gapchain=line[21]
                        if not ter:    #write TER at start of gap streak if there isnt one already
                            f.write(terline+'\n')
                            ter=True
                        if gapres==-1: #alternate gapres so each gap comes out as a different res when renaming res's
                            gapres=-2
                        elif gapres==-2:
                            gapres=-1
                        for gapatom in atomlist: 
                            f.write(gapline_func(gapchain,gapres,gapatom)+'\n')
                        col_index=col_index+1
                        gapped=True
                    else:
                        if gapped: #write ter at the end of gap streak
                            f.write(terline+'\n')
                            ter=True
                        #write line
                        prev_chain=line[21]
                        f.write(line)
                        prev_res=int(line[22:26])
                        gapped=False
                        col_index=col_index+1
                        ter=False
                        break
            elif line.startswith('TER'):
                if not ter:
                    f.write(terline+'\n')
                    ter=True
            else:
                print("potentially strange line:",line)
                f.write('\n')
                continue
        #write term gaps
        gapchain=prev_chain
        for n in range(col_index, alignment_length+1):
            if record.seq[n-1:n] == '-':
                if not ter:    #write TER at start of gap streak if there isnt one already
                    f.write(terline+'\n')
                    ter=True
                if gapres==-1: #alternate gapres so each gap comes out as a different res when renaming res's
                    gapres=-2
                elif gapres==-2:
                    gapres=-1
                for gapatom in atomlist: 
                    f.write(gapline_func(gapchain,gapres,gapatom)+'\n')
                col_index=col_index+1
                gapped=True
            else:
                emsg = 'ERROR: there are non-gap columns after the last line in pdb {} \n'.format(pdbid)
                sys.stderr.write(emsg)
                sys.exit(1)

        if gapped: #write ter at the end of gap streak
            if not ter:
                f.write(terline+'\n')
                ter=True
    else:
        continue

f.close()    
pdbfh.close()
