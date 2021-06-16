#!/usr/bin/env python3
import os
import sys
import re
from Bio import AlignIO

af = sys.argv[1] #"../mcm_test/pdb_bank/aligned/alignment_files/align_output.fasta" #sys.argv[1]
#args=["902:1101","1735:1956","2608:2703","3495:3633","4189:4289","4891:5073"]


#translate numbers
i_start_list=list()
i_stop_list=list()
align = AlignIO.read(af, "fasta") 
pdbid=os.environ["id"]
pdb=pdbid + ".pdb"
sequence=""
res_tick=0
tick=0
psuedo_n=0
translated=list()

#MIGHT NEED TO FIX TRANSLATED SECTION BELOW. IS SLIGHTLY DIFFERENT THAN CLEAN_PDB_WITHOUTGAPS.PY SECTION
for record in align:
    if record.id == pdbid:
        #get i set from residue_range
        for i in sys.argv[2:]: 
            i_split=re.split(':', i)
            i_start=int(i_split[0])-1 #-1 because splice # takes places after character #
            i_stop=int(i_split[1])
            i_start_list.append(i_start+1)
            i_stop_list.append(i_stop)
        #translate alignment column # to pdb res #
        for n in range(1, i_stop_list[-1]+1):
            subtractor=0 #the amount of gaps in the last encountered stretch of gaps
            #count residues 
            if record.seq[n-1:n] == '-':
                if n == i_start_list[tick]:
                    res_tick=res_tick+1
                else:
                    pass
            else:
                res_tick=res_tick+1
            #put important residues in pdb
            if n == i_start_list[tick]:  
                for x in range(res_tick, (res_tick+
                                (1+i_stop_list[tick]-i_start_list[tick]))):
                    if psuedo_n < n:
                        psuedo_n=n
                    if record.seq[psuedo_n-1:psuedo_n] == '-':
                        subtractor=subtractor+1
                        translated.append("GAP")
                        #pass
                    else:
                        if subtractor>=x:
                            print("WARNING: the subtractor >= x. this is extremely odd and likely indicative of an issue")
                            res=x
                        else:
                            res=x-subtractor
                        translated.append(res)
                    psuedo_n=psuedo_n+1
                tick=tick+1
                if tick+1 > len(i_start_list):
                    break
#write new pdb         
if not translated:
    print(pdbid," does not contain any residues in the specified range. Omitting from rest of process")
    os.remove(pdb)
    sys.exit()
newpdb=pdbid + "-cut" + ".pdb"
f = open(newpdb, "w")
pdbfh = open(pdb, 'r')
prev_t=translated[0]
#for x in range(len(translated)): #find first non-GAP t
#    if prev_t == 'GAP':
#        prev_t=translated[x]
#    if prev_t != 'GAP':
#        break
ter=False
prev_t_index=0
terline='TER'+' '*77
gapchain=1
gapres=-1
prev_chain=None
#gapline='ATOM      0  X   GAP 0   0       nan     nan     nan    nan   nan'
def gapline_func(gapchain,gapres):
    gapline='ATOM      0  X   GAP '+str(gapchain)+str(gapres).rjust(4)+'       nan     nan     nan    nan   nan'
    return gapline
gapped=False
for line in pdbfh:
    if line.startswith('ATOM'):
        for t in translated:
            if int(line[22:26]) == t:
                if prev_chain == None: #name chain for getchain 
                    gapchain=line[21]
                elif prev_chain != line[21]:
                    gapchain=prev_chain
                elif prev_chain == line[21]:
                    gapchain=line[21]
                if prev_t != t:
                    if prev_t != 'GAP':
                        prev_t_index=prev_t_index+1
                    gapped=False
                    for x in range(len(translated)):
                        if 'GAP' == translated[prev_t_index]:
                            if not ter:    #prevents double TER
                                f.write(terline+'\n')
                                ter=True
                            if gapres==-1: #alternate gapres so each gap comes out as a different res when renaming res's
                                gapres=-2
                            elif gapres==-2:
                                gapres=-1
                            f.write(gapline_func(gapchain,gapres)+'\n')
                            prev_t_index=prev_t_index+1
                            gapped=True
                    if gapped: #write ter at the end of gap streak
                            f.write(terline+'\n')
                            ter=True
                    if prev_t != 'GAP':
                        if t != prev_t+1:
                                if not ter:    #prevents double TER
                                    f.write(terline+'\n')
                                    ter=True
                                    realter=False
                prev_t=t
                prev_chain=line[21]
                f.write(line)
                ter=False
                realter=False
            else:
                continue
    elif line.startswith('TER'):
        if not ter:
            f.write(terline+'\n')
            realter=True
            ter=True
    else:
        print("potentially strange line:",line)
        f.write('\n')

#write terminal gap streak
if prev_t != 'GAP':
    prev_t_index=prev_t_index+1
for x in range(len(translated)):
    try:
        if 'GAP' == translated[prev_t_index]:
            if not ter:    #prevents double TER
                f.write(terline+'\n')
                ter=True
            if gapres==-1: #alternate gapres so each gap comes out as a different res when renaming res's
                gapres=-2
            elif gapres==-2:
                gapres=-1
            f.write(gapline_func(gapchain,gapres)+'\n')
            prev_t_index=prev_t_index+1
            gapped=True
    except:
        break
if gapped: #write ter at the end of gap streak
        f.write(terline+'\n')
        ter=True
        
f.close()    
pdbfh.close()

