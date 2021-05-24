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


#make unmade directories and csvs
mkdir -p $label/{pdb_bank/{crude,aligned,cut_and_cleaned},blast_files,metadata}
cd $label
fasta="../$fasta"
printf 'select,pdb,chain,chains,description,species,residstart,end,ligand,mutation,resolution,methods,complex,DNA,nucleotide,nt type,intpartner,CA,remark \n' > metadata/${label}_mastertable.csv
echo "made $label/metadata/${label}_mastertable.csv"
for m in ${compare_mer[@]}; do
	mkdir -p pdb_bank/cut_and_cleaned/$m-mer
	printf 'select,pdb,chain,species,%s-mers,residstart,end,DNA,nucleotide,dist,chainz,number,int,gap \n' "$m" > metadata/${label}_svd_${m}-mer.csv
	echo "made $label/metadata/${label}_svd_${m}-mer.csv"
done

#write log
log="$label".log
printf '\n%s' "Step 3 for $label - " 2>&1 | tee -a $log
date 2>&1 | tee -a $log
echo "Step 3 Variables:"
printf '\tComparing oligomerization states:'; print_array "${compare_mer[@]}"
printf '\tResidue Range:'; print_array "${residue_range[@]}"

#cut pdbs according to alignment 
header trimming and cleaning pdbs
cp pdb_bank/aligned/standard_out/*drd.pdb pdb_bank/cut_and_cleaned
cd pdb_bank/cut_and_cleaned
#cp ../aligned/pdb_list_align_input.dat pdb_list_clean.dat
u_align=$(wc -l ../aligned/pdb_list_align_input.dat | awk '{ print $1 }')
rm -f ../cut_and_cleaned/pdb_list_clean.dat
for row in $(seq 1 "$u_align") ; do
	id=$(awk "FNR == $row { print \$1 }" ../aligned/pdb_list_align_input.dat)
	echo "$id" 2>&1 | tee -a ../../$log
	#rename chains
	echo "	renaming chains to letters" 2>&1 | tee -a ../../$log
	pdb_chainbows $id.pdb > $id-chains.pdb
	rm $id.pdb
	mv $id-chains.pdb $id.pdb
	#trim pdb
	echo "	trimming" 2>&1 | tee -a ../../$log
	python3 ../../../code/clean_pdb.py ${residue_range[@]} 2>&1 | tee -a ../../$log
	if [ ! -f "$id.pdb" ]; then
  	      	echo "moving on to the next pdb..." 2>&1 | tee -a ../../$log
		continue
	fi
	rm $id.pdb
	mv $id-cut.pdb $id.pdb
	#renumber atoms
	echo "	renumbering atoms" 2>&1 | tee -a ../../$log
	pdb_reatom -1 $id.pdb > $id-tmp.pdb
	rm $id.pdb
	mv $id-tmp.pdb $id.pdb
	#renumber residues
	echo "	renumbering residues" 2>&1 | tee -a ../../$log
	pdb_reres -1 $id.pdb > $id-tmp.pdb
	rm $id.pdb
	mv $id-tmp.pdb $id.pdb
	#reorder chains starting from A
		#placeholder in case i need to do this. doesnt work 
#	for alpha in {A,B,C,D,E,F}; do	
#		pdb_reres rplchain -A:$alpha $id.pdb > $id-tmp.pdb
#		rm $id.pdb
#		mv $id-tmp.pdb $id.pdb
#
#	done
	echo "$id" >> ../cut_and_cleaned/pdb_list_clean.dat
done


cd ../..
footer trimmed and cleaned pdbs

#rename chains based on homology # in step 2 now
#header standardizing chains ${pdb_ref_chains[@]} according to $pdb_ref_id 
#cd pdb_bank/standardized_chains
#python3 ../../../code/align_chains.py ${pdb_ref_chains[@]}
#cd ../..
#footer standardized chains ${pdb_ref_chains[@]} according to $pdb_ref_id

#sep pdbs into mers 
header separating pdbs into n-mers
cd pdb_bank/cut_and_cleaned
u_clean=$(wc -l ../cut_and_cleaned/pdb_list_clean.dat | awk '{ print $1 }')
for row2 in $(seq 1 "$u_clean") ; do
	printf '\n'

	id=$(awk "FNR == $row2 { print \$1 }" ../cut_and_cleaned/pdb_list_clean.dat)
	if [ ! -f "$id.pdb" ]; then
		echo "$id.pdb was not found. moving on to the next one" 2>&1 | tee -a ../../$log
		continue
	fi
	
	
	#get total chains
	chains=$(pdb_wc $id.pdb | awk "FNR == 2 { print \$3 }")
	chn_ar=() #reset
	for chn_nbr in $(seq 1 "$chains") ; do
		chn_ltr=`python3 -c "print(chr(ord('A')+$(($chn_nbr-1))))"`
		chn_ar[$chn_nbr]=$chn_ltr
	done
	allchains=$(echo ${chn_ar[@]} | tr -d ' ')
	printf "yes,$id,$allchains \n" >> ../../metadata/${label}_mastertable.csv
	
	#split for each oligo state
	if [[ " ${compare_mer[@]} " =~ " 1 " ]]; then
		mer=1
		cp $id.pdb $mer-mer
		cd $mer-mer
		echo "separating $id into $mer-mers..." 2>&1 | tee -a ../../../$log
		if [ $chains -lt $mer ] ; then
			echo "$id has less than $mer chains (it has $chains) so it cant be used for $mer-mer comparisons." 2>&1 | tee -a ../../../$log
			rm $id.pdb
			cd ..
			continue
		fi
		echo "from chains $allchains, wrote:" 2>&1 | tee -a ../../../$log
		for chn in $(seq 1 "$chains"); do
			set +e
			pdb_selchain -${chn_ar[$chn]} $id.pdb > ${id}_${chn_ar[$chn]}.pdb
			checkchains ${id}_${chn_ar[$chn]}.pdb $mer
			echo "    ${id}_${chn_ar[$chn]}.pdb" 2>&1 | tee -a ../../../$log
			printf "yes,${id}_${chn_ar[$chn]},${chn_ar[$chn]},,${chn_ar[$chn]} \n" >> ../../../metadata/${label}_svd_$mer-mer.csv
			set -e
		done
		rm $id.pdb
		cd ..
	fi
	
	splitpdb 2 ${dimer_pattern[@]}
	splitpdb 4 ${tetramer_pattern[@]}
	splitpdb 6 ${hexamer_pattern[@]}

done
cd ../..
footer separated pdbs into n-mers


(printf '%0.s^v' {1..50}; printf '\n')  >> $log
(printf '%0.sv^' {1..50}; printf '\n')  >> $log
cd ..
