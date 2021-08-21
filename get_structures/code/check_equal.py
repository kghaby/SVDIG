#!/usr/bin/env python3
import os
import sys
#sys.argv[1]
#Checking if any residue contains more atoms than specified in master.sh.
#print('\t',"Checking if any residue contains more atoms than specified in master.sh ({})".format(os.environ["atomcount"]))
print('\tchecking that atoms are equal')
pdbid = os.environ["newid"] 
pdbf = pdbid+".pdb"
pdbfh = open(pdbf, 'r')
prev_res=None
count=1
overcount=0
atomcount=int(os.environ["atomcount"])
badres=set()
extraatoms=list()
badresindex=0
        
for line in pdbfh:
    if line.startswith('ATOM'):
        if line[22:26] != prev_res:
            count=1
        if line[22:26] == prev_res:
            count=count+1
            if count > atomcount: #if the residue has more lines than there are atoms in the master.sh "atoms" variable
                badres.add(int(line[22:26]))
                if len(extraatoms) < len(badres):
                    extraatoms.append(1)
                else:
                    extraatoms[badresindex]=extraatoms[badresindex]+1
        prev_res=line[22:26]
        badresindex=len(badres)-1 #-1 bc 0 is first position in list, not 1 

if len(badres) == 0:
    sys.exit()

badres_l=list(badres)

if len(badres) > 0:
    print('\t\t',"WARNING: Found {} residues with more atoms than what was specified in master.sh ({})".format(len(badres),atomcount))
for x in range(len(badres)):
    print('\t\t\t',"Residue {} had {} extra atoms".format(badres_l[x],extraatoms[x]))

#reopen
pdbfh.close()
pdbfh = open(pdbf, 'r')

pdbfh_new=open("{}_fixed.pdb".format(pdbid), 'w')

#fix 
#detect if ares and bres 

print('\t\t',"Attempting to identify duplicates styled as AXXX and BXXX")

#HAVE TO IDENTIFY THAT NEXT LINE IS DUP BC SOME ATOM NAMES MAY BLEED INTO 16:17. CANT JUST LOOK AT THAT POSITION

for line in pdbfh:
    if line.startswith('ATOM'):
        if int(line[22:26]) in badres_l:
            if line[16:17] != " ":
                if line[16:17] != "A":
                    #print(line[16:17])
                    print('\t\t\t',"DELETING LINE:")
                    print(line)
                    continue
                else:
                    pdbfh_new.write(line)
            else: 
                pdbfh_new.write(line)
        else: 
            pdbfh_new.write(line)
                        
    else:
        pdbfh_new.write(line)
        
