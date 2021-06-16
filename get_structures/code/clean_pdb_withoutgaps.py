#!/usr/bin/env python3
import os
import sys
import re
from Bio import AlignIO

af = sys.argv[1]
align = AlignIO.read(af, "fasta")

#translate numbers
i_start_list=list()
i_stop_list=list()
pdbid=os.environ["id"]
pdb=pdbid + ".pdb"
sequence=""
res_tick=0
tick=0
psuedo_n=0
translated=list()
for record in align:
    if record.id == pdbid:
        #get i set
        for i in sys.argv[2:]:
            i_split=re.split(':', i)
            i_start=int(i_split[0])
            i_stop=int(i_split[1])
            i_start_list.append(i_start)
            i_stop_list.append(i_stop)

        #translate alignment column # to pdb res #
        for n in range(1, i_stop_list[-1]+1):
            subtractor=0
            #count residues
            if record.seq[n-1] == '-':
                pass
            else:
                res_tick=res_tick+1
            #put important residues in pdb
            if n == i_start_list[tick]:
                for x in range(res_tick, (res_tick+(i_stop_list[tick]-(i_start_list[tick]))+1)):
                    if psuedo_n < n:
                        psuedo_n=n
                    if record.seq[psuedo_n-1] == '-':
                        if n != psuedo_n:
                            subtractor=subtractor+1
                        pass
                    else:
                        if subtractor>=x:
                            res=x
                            print("WARNING: the subtractor >= x. this is extremely odd and likely indicative of an issue")
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
prev_writ_t=translated[0]-1
prev_t=translated[0]
ter=False
terline='TER'+' '*77
for line in pdbfh:
    if line.startswith('ATOM'):
        for t in translated:
            if int(line[22:26]) == t:
                if prev_t != t:
                    if t != prev_t+1:
                        if not ter:    #prevents double TER
                            f.write(terline+'\n')
                            ter=True
                prev_t=t
                f.write(line)
                ter=False
            else:
                continue
    elif line.startswith('TER'):
        if not ter:
            f.write(terline+'\n')
            ter=True
    else:
        print("potentially strange line:",line)
        f.write('\n')
f.close()    
pdbfh.close()
