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
printf '\n%s' "Step 3 for $label - " 2>&1 | tee -a $log
date 2>&1 | tee -a $log
echo "Step 3 Variables: " 2>&1 | tee -a $log
printf '\tReference pdb for chain standardization: %s\n' "$pdb_ref_id" 2>&1 | tee -a $log


header 'standardizing chains and aligning'
#get new number of unique pdbs
cd pdb_bank/aligned
u_std=$(wc -l pdb_list_standard_input.dat | awk '{ print $1 }')
#renumber residues in pdbs
for row in $(seq 1 "$u_std") ; do
    id=$(awk "FNR == $row { print \$1 }" pdb_list_standard_input.dat)
    mv -f $id.pdb $id.tmp.pdb
    python3 ../../../code/renumber_res.py 1 $id.tmp.pdb > $id.pdb
    rm -f $id.tmp.pdb
done

#rename chains based on homology
mkdir -p {standard_in,standard_out}
printf '\n'
echo "standardizing chains according to $pdb_ref_id" 2>&1 | tee -a ../../../$log
cd standard_in
mv ../pdb_list_standard_input.dat .
mv ../*.pdb .
python3 ../../../../code/align_chains_v2.py 2>&1 | tee -a ../../../$log
mv *_stndrd.pdb ../standard_out
cd ../standard_out
#make list of standardized pdbs
rm -f ../pdb_list_align_input.dat
for f in *_stndrd.pdb; do
    echo "${f%.pdb}" >> ../pdb_list_align_input.dat
    ter_first "${f}"
    ter_last "${f}"
done
echo "standardized chains according to $pdb_ref_id"
printf '\n'

#get new number of unique pdbs
u_align=$(wc -l ../pdb_list_align_input.dat | awk '{ print $1 }')
#renumber residues AND atoms
    #these still include unwanted residues, so it will be renumbered again. these numbers are a reference for alignment
for row in $(seq 1 "$u_align") ; do
    id=$(awk "FNR == $row { print \$1 }" ../pdb_list_align_input.dat)
    mv -f $id.pdb $id.tmp.pdb
    python3 ../../../../code/renumber_res.py 1 $id.tmp.pdb > $id.tmp2.pdb
    python3 ../../../../code/renumber_atoms.py 1 $id.tmp2.pdb > $id.pdb
    rm -f $id.tmp.pdb
    rm -f $id.tmp2.pdb
done

#get seqs for each pdb and do pairwise alignment from input seq to each new seq
for row in $(seq 1 "$u_align") ; do
    id=$(awk "FNR == $row { print \$1 }" ../pdb_list_align_input.dat)
    pdb_tofasta $id.pdb > $id.fasta
    newfasta="$id.fasta"
    #python3 ../../../code/align_pairwise.py > ${id}_alignment.fasta #placeholder. doesnt work bc isnt needed atm
done
mv -f *fasta ..
cd ..
printf '\n'

#align all
set +e
echo "aligning now...will take a sec"
python3 ../../../code/align_all.py pdb_list_align_input.dat
cp ../../$fasta .
mkdir -p alignment_files
mv -f *fasta alignment_files
echo "if using mafft: open $label/pdb_bank/aligned/alignment_files/align_output.fasta in a viewer like Jalview. 
if not using mafft: input $label/pdb_bank/aligned/alignment_files/align_input.fasta. 
	make sure the alignment output is $label/pdb_bank/aligned/alignment_files/align_output.fasta and in fasta format" 2>&1 | tee -a ../../$log
rm -f pdb_list_mdl_split.dat 
cd ../..
footer aligned

(printf '%0.s^v' {1..50}; printf '\n')  >> $log
(printf '%0.sv^' {1..50}; printf '\n')  >> $log
cd ..
