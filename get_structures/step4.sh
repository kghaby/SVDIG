#!/bin/bash
set -ea; trap 'echo "Mission failed, we'"'"'ll get line $LINENO next time."' ERR
source main.sh

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


#add 1 to compare_mer if not already in there
if [[ ! " ${compare_mer[@]} " =~ " 1 " ]]; then
    compare_mer+=(1)
fi 

#make unmade directories and csvs
mkdir -p $label/{pdb_bank/{crude,aligned,cut_and_cleaned/withoutgaps},blast_files,metadata}
cd $label
cp ../${alignment_file} pdb_bank/cut_and_cleaned/withoutgaps/
af="${alignment_file##*/}"
fasta="../$fasta"
printf 'select,pdb,chain,chains,description,species,residstart,end,ligand,mutation,resolution,methods,complex,DNA,nucleotide,nt type,intpartner,CA,remark \n' > metadata/${label}_mastertable.csv
echo "made $label/metadata/${label}_mastertable.csv"
for m in ${compare_mer[@]}; do
	mkdir -p pdb_bank/cut_and_cleaned/withoutgaps/$m-mer
	printf 'select,pdb,chain,species,%s-mers,residstart,end,DNA,nucleotide,dist,chainz,number,int,gap \n' "$m" > metadata/${label}_svd_${m}-mer.csv
	echo "made $label/metadata/${label}_svd_${m}-mer.csv"
done

#write log
log="$label".log
printf '\n%s' "Step 4 for $label - " 2>&1 | tee -a $log
date 2>&1 | tee -a $log
echo "Step 4 Variables:"
printf '\tAlignment output file: ' 2>&1 | tee -a $log ; echo "${alignment_file}" 2>&1 | tee -a $log
printf '\tComparing oligomerization states: ' 2>&1 | tee -a $log ; print_array "${compare_mer[@]}" 2>&1 | tee -a $log
printf '\tResidue Range: ' 2>&1 | tee -a $log ; print_array "${residue_range[@]}" 2>&1 | tee -a $log


#cut pdbs according to alignment 
header trimming and cleaning pdbs 
cp pdb_bank/aligned/standard_out/*stndrd.pdb pdb_bank/cut_and_cleaned/withoutgaps
cd pdb_bank/cut_and_cleaned/withoutgaps
u_align=$(wc -l ../../aligned/pdb_list_align_input.dat | awk '{ print $1 }')
rm -f pdb_list_clean.dat
for row in $(seq 1 "$u_align") ; do
	id=$(awk "FNR == $row { print \$1 }" ../../aligned/pdb_list_align_input.dat)
	echo "$id" 2>&1 | tee -a ../../../$log
	#rename chains
	echo "	renaming chains to letters" 2>&1 | tee -a ../../../$log
	python3 ../../../../code/reletter_chains.py $id.pdb > $id-chains.pdb
	rm -f $id.pdb
	mv $id-chains.pdb $id.pdb
    ter_first "${id}.pdb"
    #trim pdb
	echo "	trimming" 2>&1 | tee -a ../../../$log
	python3 ../../../../code/clean_pdb_withoutgaps.py $af ${residue_range[@]} 2>&1 | tee -a ../../../$log
	if [ ! -f "$id.pdb" ]; then
  	    printf '\t' 2>&1 | tee -a ../../../$log
        echo "moving on to the next pdb..." 2>&1 | tee -a ../../../$log
		continue
	fi
	rm -f $id.pdb
	mv $id-cut.pdb $id.pdb
    ter_first "${id}.pdb"
    ter_last "${id}.pdb" 
    #renumber atoms
	echo "	renumbering atoms" 2>&1 | tee -a ../../../$log
	python3 ../../../../code/renumber_atoms.py 1 $id.pdb > $id-tmp.pdb
	rm -f $id.pdb
	mv $id-tmp.pdb $id.pdb
	#renumber residues
	echo "	renumbering residues" 2>&1 | tee -a ../../../$log
	python3 ../../../../code/renumber_res.py 1 $id.pdb > $id-tmp.pdb
	rm -f $id.pdb
	mv $id-tmp.pdb $id.pdb
	echo "$id" >> pdb_list_clean.dat
done
cd ../../../
footer trimmed and cleaned pdbs


#sep pdbs into mers 
header separating pdbs into n-mers
cd pdb_bank/cut_and_cleaned/withoutgaps
for mer in ${compare_mer[@]}; do
    rm -f ${mer}-mer/pdb_list_${mer}-mer.dat
done

u_clean=$(wc -l pdb_list_clean.dat | awk '{ print $1 }')
for row2 in $(seq 1 "$u_clean") ; do
	printf '\n'

	id=$(awk "FNR == $row2 { print \$1 }" pdb_list_clean.dat)
	if [ ! -f "$id.pdb" ]; then
		echo "$id.pdb was not found. moving on to the next one" 2>&1 | tee -a ../../../$log
		continue
	fi
	
	#get total chains
	chains=$(pdb_wc $id.pdb | awk "FNR == 2 { print \$3 }")
	chn_ar=() #reset

    #this chunk just prints a seq of letters based off of the number of chains. it doesnt account for a non-A start. thus its bad
        #keeping here bc it may be useful one day 
#	for chn_nbr in $(seq 1 "$chains") ; do
#		chn_ltr=`python3 -c "print(chr(ord('A')+$(($chn_nbr-1))))"`
#		chn_ar[$chn_nbr]=$chn_ltr
#	done
    
    chn_ar=($(python3 ../../../../code/identify_chains.py | tr -d "[],'"))
	allchains=$(echo ${chn_ar[@]} | tr -d ' ')
	printf "yes,$id,$allchains \n" >> ../../../metadata/${label}_mastertable.csv
	
	#split for each oligo state 
    #1-mer if-then block
	if [[ " ${compare_mer[@]} " =~ " 1 " ]]; then
		mer=1
		cp $id.pdb $mer-mer
		cd $mer-mer
		echo "separating $id into $mer-mers..." 2>&1 | tee -a ../../../../$log
		if [ $chains -lt $mer ] ; then
			echo "$id has less than $mer chains (it has $chains) so it cant be used for $mer-mer comparisons." 2>&1 | tee -a ../../../../$log
			rm -f $id.pdb
			cd ..
			continue
		fi
		echo "from chains $allchains, wrote:" 2>&1 | tee -a ../../../../$log
		for chn in $(seq 0 "$((${chains}-1))"); do
			set +e
			python3 ../../../../../code/take_chains.py -${chn_ar[$chn]} $id.pdb > ${id}_${chn_ar[$chn]}.pdb
            ter_first "${id}_${chn_ar[$chn]}.pdb"
            ter_last "${id}_${chn_ar[$chn]}.pdb"
            checkchains ${id}_${chn_ar[$chn]}.pdb $mer
            echo "${id}_${chn_ar[$chn]}" >> pdb_list_${mer}-mer.dat
			echo "    ${id}_${chn_ar[$chn]}.pdb" 2>&1 | tee -a ../../../../$log
			printf "yes,${id}_${chn_ar[$chn]},${chn_ar[$chn]},,${chn_ar[$chn]} \n" >> ../../../../metadata/${label}_svd_$mer-mer.csv
			set -e
		done
		rm -f $id.pdb
		cd ..
    fi
	
    splitpdb 2 ${dimer_pattern[@]}
    splitpdb 4 ${tetramer_pattern[@]}
    splitpdb 6 ${hexamer_pattern[@]}

done
cd ../../../
footer separated pdbs into n-mers


(printf '%0.s^v' {1..50}; printf '\n')  >> $log
(printf '%0.sv^' {1..50}; printf '\n')  >> $log
cd ..
