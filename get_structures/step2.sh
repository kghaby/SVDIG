#!/bin/bash
set -ea; trap 'echo "Mission failed, we'"'"'ll get line $LINENO next time."' ERR
source master.sh

#MAIN
#==================================================

#check if label has been used
if [ ! -d "$label" ]; then
	echo "The label $label has not been used. You should start with step1.sh"
	select yn in "Use existing directory" "Exit script"; do
		case $yn in
			"Use existing directory" ) echo "Moving on..."; break;;
			"Exit script" ) echo "Change label at the top of master.sh"; exit;;
		esac
	done
fi

if [ -d "$label" ]; then
    echo "To confirm, use the directory labeled $label?"  
    select yn in "Use existing directory" "Exit script"; do
        case $yn in
            "Use existing directory" ) echo "Moving on."; break;;
            "Exit script" ) echo "Change label at the top of this file"; exit;;
        esac
    done
fi

#make unmade directories
mkdir -p $label/{pdb_bank/{crude,aligned,cut_and_cleaned},blast_files,metadata}
cd $label
fasta="../$fasta"

#write log
log="$label".log
printf '\n%s' "Step 2 for $label - " 2>&1 | tee -a $log
date 2>&1 | tee -a $log
echo "Step 2 Variables: " 2>&1 | tee -a $log
printf '\tE thresh for rest of process: %s\n' "$newethresh" 2>&1 | tee -a $log

#get crude pdbs
header downloading pdbs
python3 ../code/id_pdbs.py > pdb_bank/crude/pdb_list.dat 
uniques=$(wc -l pdb_bank/crude/pdb_list.dat | awk '{ print $1 }')
cd pdb_bank/crude
for row in $(seq 1 "$uniques") ; do
	printf '\n'
	id=$(awk "FNR == $row { print \$1 }" ../crude/pdb_list.dat)
	echo "looking for $id.pdb" 2>&1 | tee -a ../../$log
	python3 ../../../code/fetch_pdbs.py 2>&1 | tee -a ../../$log 
	if [ -f "$id.pdb" ]; then
		echo "got it. heres a summary of $id's contents:" 2>&1 | tee -a ../../$log
		pdb_wc $id.pdb 2>&1 | tee -a ../../$log
	fi
done
cd ../..
footer downloaded pdbs

#formatting stuff
header 'parsing pdbs'
cd pdb_bank/aligned
cp ../crude/*.pdb .
reffasta="../../$fasta" 
#split models into own pdbs
rm -f pdb_list_mdl_split.dat
for row in $(seq 1 "$uniques") ; do
	id=$(awk "FNR == $row { print \$1 }" ../crude/pdb_list.dat)
	skip_if_not_found "$id.pdb"
	models=$(pdb_wc $id.pdb | awk "FNR == 1 { print \$3 }")
	if [ $models -gt 1 ]; then
		echo "PDB $id has $models models. splitting models into their own pdbs" 2>&1 | tee -a ../../$log
		pdb_splitmodel $id.pdb
		for mdlnbr in $(seq 1 "$models"); do
			echo "${id}_$mdlnbr" >> pdb_list_mdl_split.dat
			rm -f $id.pdb
		done
		continue
	fi
	echo "$id" >> pdb_list_mdl_split.dat
done
u_mdl_split=$(wc -l ../aligned/pdb_list_mdl_split.dat | awk '{ print $1 }')
#take atom lines
for row in $(seq 1 "$u_mdl_split") ; do
	id=$(awk "FNR == $row { print \$1 }" ../aligned/pdb_list_mdl_split.dat)
	echo "grabbing atom lines in $id" 2>&1 | tee -a ../../$log
	python3 ../../../code/atoms.py $id.pdb > $id-atoms.pdb
	rm -f $id.pdb
done
#split numbered and lettered chains
rm -f pdb_list_standard_input.dat
for row in $(seq 1 "$u_mdl_split") ; do
	id=$(awk "FNR == $row { print \$1 }" ../aligned/pdb_list_mdl_split.dat)
	chnfile="$id-atoms.pdb"
	chnseq=$(python3 ../../../code/id_chains.py)
	if [[ $chnseq == *[[:digit:]]* ]] && [[ $chnseq == *[[:alpha:]]* ]]; then
		echo "PDB $id contains lettered and numbered chains. splitting them into their own pdbs, ${id}_ltr.pdb and ${id}_nbr.pdb"
		chnseqltr=${chnseq//[[:digit:]]/}
		chnseqltr_comma=$(echo $chnseqltr | sed -r 's/(.{1})/\1,/g;s/,$//') #put commas in between chains
		python3 ../../../code/take_chains.py -$chnseqltr_comma $id-atoms.pdb > ${id}_ltr.pdb
        ter_first "${id}_ltr.pdb"
        ter_last "${id}_ltr.pdb"
		chnseqnbr=${chnseq//[[:alpha:]]/}
		chnseqnbr_comma=$(echo $chnseqnbr | sed -r 's/(.{1})/\1,/g;s/,$//') #put commas in between chains
		python3 ../../../code/take_chains.py -$chnseqnbr_comma $id-atoms.pdb > ${id}_nbr.pdb
        ter_first "${id}_nbr.pdb"
        ter_last "${id}_nbr.pdb"
        echo "${id}_ltr" >> pdb_list_standard_input.dat
		echo "${id}_nbr" >> pdb_list_standard_input.dat
		rm $id-atoms.pdb
		continue	 		
	fi 
	mv $id-atoms.pdb $id.pdb
	echo "$id" >> pdb_list_standard_input.dat
done
cd ../..
footer parsed pdbs

echo "pdb ids are in $label/pdb_bank/aligned/pdb_list_standard_input.dat so you may use this file to choose a chain reference pdb for step 3 if you dont already know it" 2>&1 | tee -a $log

(printf '%0.s^v' {1..50}; printf '\n')  >> $log
(printf '%0.sv^' {1..50}; printf '\n')  >> $log
cd ..
