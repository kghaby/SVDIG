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


#make unmade directories
mkdir -p $label/{pdb_bank/{crude,aligned,cut_and_cleaned},blast_files,metadata}
cd $label
fasta="../$fasta"

#write log
log="$label".log
printf '\n%s' "Step 2 for $label - " 2>&1 | tee -a $log
date 2>&1 | tee -a $log
echo "Step 2 Variables:"
printf '\tE thresh for rest of process: %s\n' "$newethresh"

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

#alignment 
header aligning
cd pdb_bank/aligned
cp ../crude/*.pdb .
reffasta="../../$fasta" 
#split models into own pdbs
rm -f ../aligned/pdb_list_mdl_split.dat
for row in $(seq 1 "$uniques") ; do
	id=$(awk "FNR == $row { print \$1 }" ../crude/pdb_list.dat)
	skip_if_not_found "$id.pdb"
	models=$(pdb_wc $id.pdb | awk "FNR == 1 { print \$3 }")
	if [ $models -gt 1 ]; then
		echo "PDB $id has $models models. splitting models into their own pdbs" 2>&1 | tee -a ../../$log
		pdb_splitmodel $id.pdb
		for mdlnbr in $(seq 1 "$models"); do
			echo "${id}_$mdlnbr" >> ../aligned/pdb_list_mdl_split.dat
			rm $id.pdb
		done
		continue
	fi
	echo "$id" >> ../aligned/pdb_list_mdl_split.dat
done
u_mdl_split=$(wc -l ../aligned/pdb_list_mdl_split.dat | awk '{ print $1 }')
#take atom lines
for row in $(seq 1 "$u_mdl_split") ; do
	id=$(awk "FNR == $row { print \$1 }" ../aligned/pdb_list_mdl_split.dat)
	echo "grabbing atom lines in $id" 2>&1 | tee -a ../../$log
	python3 ../../../code/atoms.py $id.pdb > $id-atoms.pdb
	rm $id.pdb
done
#split numbered and lettered chains
rm -f ../aligned/pdb_list_align_input.dat
for row in $(seq 1 "$u_mdl_split") ; do
	id=$(awk "FNR == $row { print \$1 }" ../aligned/pdb_list_mdl_split.dat)
	chnfile="$id-atoms.pdb"
	chnseq=$(python3 ../../../code/id_chains.py)
	if [[ $chnseq == *[[:digit:]]* ]] && [[ $chnseq == *[[:alpha:]]* ]]; then
		echo "PDB $id contains lettered and numbered chains. splitting them into their own pdbs, ${id}_ltr.pdb and ${id}_nbr.pdb"
		chnseqltr=${chnseq//[[:digit:]]/}
		chnseqltr_comma=$(echo $chnseqltr | sed -r 's/(.{1})/\1,/g;s/,$//') #put commas in between chains
		pdb_selchain -$chnseqltr_comma $id-atoms.pdb > ${id}_ltr.pdb
		chnseqnbr=${chnseq//[[:alpha:]]/}
		chnseqnbr_comma=$(echo $chnseqnbr | sed -r 's/(.{1})/\1,/g;s/,$//') #put commas in between chains
		pdb_selchain -$chnseqnbr_comma $id-atoms.pdb > ${id}_nbr.pdb
		echo "${id}_ltr" >> ../aligned/pdb_list_align_input.dat
		echo "${id}_nbr" >> ../aligned/pdb_list_align_input.dat
		rm $id-atoms.pdb
		continue	 		
	fi 
	mv $id-atoms.pdb $id.pdb
	echo "$id" >> ../aligned/pdb_list_align_input.dat
done

#get new number of unique pdbs
u_align=$(wc -l ../aligned/pdb_list_align_input.dat | awk '{ print $1 }')

#renumber residues in pdbs 
	#these still include unwanted residues, so it will be renumbered again. these numbers are a reference for alignment
for row in $(seq 1 "$u_align") ; do
	id=$(awk "FNR == $row { print \$1 }" ../aligned/pdb_list_align_input.dat)
	mv -f $id.pdb $id.tmp.pdb
	pdb_reres -1 $id.tmp.pdb > $id.pdb
	rm $id.tmp.pdb  
done

#get seqs for each pdb and do pairwise alignment from input seq to each new seq
for row in $(seq 1 "$u_align") ; do
	id=$(awk "FNR == $row { print \$1 }" ../aligned/pdb_list_align_input.dat)
	pdb_tofasta $id.pdb > $id.fasta
	newfasta="$id.fasta"
	#python3 ../../../code/align_pairwise.py > ${id}_alignment.fasta #placeholder. doesnt work bc isnt needed atm 
done

#align all
echo "aligning now...will take a sec"
python3 ../../../code/align_all.py
cp ../../$fasta .
mkdir -p alignment_files
mv -f *fasta alignment_files
echo "if using mafft: open $label/pdb_bank/aligned/alignment_files/align_output.fasta in a viewer like Jalview. 
if not using mafft: input $label/pdb_bank/aligned/alignment_files/align_input.fasta. 
	make sure the alignment output is $label/pdb_bank/aligned/alignment_files/align_output.fasta and in fasta format" 2>&1 | tee -a ../../$log
rm ../aligned/pdb_list_mdl_split.dat 
cd ../..
footer aligned

(printf '%0.s^v' {1..50}; printf '\n')  >> $log
(printf '%0.sv^' {1..50}; printf '\n')  >> $log
cd ..
