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
    echo "To confirm, use the directory labeled $label? CSVs from step 4 will not be touched by the way."
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

#make unmade directories
mkdir -p ${label}/{pdb_bank/equal_atoms,metadata}
cd $label
for m in ${compare_mer[@]}; do
	mkdir -p pdb_bank/equal_atoms/$m-mer
done

#write log
log="$label".log
printf '\n%s' "Step 5 for $label - " 2>&1 | tee -a $log
date 2>&1 | tee -a $log
echo "Step 5 Variables: " 2>&1 | tee -a $log
printf '\tComparing oligomerization states: ' 2>&1 | tee -a $log ; print_array "${compare_mer[@]}" 2>&1 | tee -a $log
printf '\tSelecting atoms: ' 2>&1 | tee -a $log ; print_array "${atoms[@]}" 2>&1 | tee -a $log


#get builder pdbs
header writing builder pdbs with equal atoms
m=1
pdb_list_s4_input="pdb_bank/cut_and_cleaned/withoutgaps/${m}-mer/pdb_list_${m}-mer.dat"
cp $pdb_list_s4_input pdb_bank/equal_atoms/${m}-mer/
cd pdb_bank/equal_atoms/${m}-mer/
printf 'pdb-id\tatoms\tresidues\tchains\n' > ${m}-mer_info.txt
#get seqs
echo "getting sequences of ${m}-mers" 2>&1 | tee -a ../../../$log
mkdir -p alignment_files
u_mer=$(wc -l pdb_list_${m}-mer.dat | awk '{ print $1 }')
for row in $(seq 1 "$u_mer") ; do
    id=$(awk "FNR == $row { print \$1 }" pdb_list_${m}-mer.dat)
    cp ../../cut_and_cleaned/withoutgaps/${m}-mer/$id.pdb .
    pdb_tofasta $id.pdb > $id.fasta
done
#align
echo "aligning ${m}-mers" 2>&1 | tee -a ../../../$log
python3 ../../../../code/align_all.py pdb_list_${m}-mer.dat
mv *.fasta alignment_files
#equalize_atoms
for row in $(seq 1 "$u_mer") ; do
    id=$(awk "FNR == $row { print \$1 }" pdb_list_${m}-mer.dat)
    echo "${id}.pdb" 2>&1 | tee -a ../../../$log
    printf '\t' 2>&1 | tee -a ../../../$log
    echo "equalizing atoms" 2>&1 | tee -a ../../../$log
    python3 ../../../../code/equalize_atoms.py alignment_files/align_output.fasta ${atoms[@]} 2>&1 | tee -a ../../../$log
    rm -f ${id}.pdb

    ter_first "${id}_withgaps.pdb"
    ter_last "${id}_withgaps.pdb"

    printf '\t' 2>&1 | tee -a ../../../$log
    echo "reordering atoms, residues, and chains" 2>&1 | tee -a ../../../$log
    #renumber atoms
    python3 ../../../../code/renumber_atoms.py 1 ${id}_withgaps.pdb > ${id}_a.pdb 
    rm -f ${id}_withgaps.pdb
    #renumber residues
    python3 ../../../../code/renumber_res.py 1 ${id}_a.pdb > ${id}_withgaps.pdb 
    rm -f ${id}_a.pdb
    
    #check atoms
    #echo "checking atoms" 2>&1 | tee -a ../../../$log
    newid="${id}_withgaps"
    atomcount=${#atoms[@]} 
    python3 ../../../../code/check_equal.py 2>&1 | tee -a ../../../$log 
    if [ -f "${newid}_fixed.pdb" ]; then
        python3 ../../../../code/renumber_atoms.py 1 ${newid}_fixed.pdb > ${newid}.pdb  
        rm -f ${newid}_fixed.pdb 
    fi
    
    #write info
    a=$(pdb_wc -a ${id}_withgaps.pdb | awk "FNR == 1 { print \$3 }")
    r=$(pdb_wc -r ${id}_withgaps.pdb | awk "FNR == 1 { print \$3 }")
    c=$(pdb_wc -c ${id}_withgaps.pdb | awk "FNR == 1 { print \$3 }")
    printf '%s\t%s\t%s\t%s\n' "${id}_withgaps.pdb" "$a" "$r" "$c" >> ${m}-mer_info.txt
   
    printf '\n' 2>&1 | tee -a ../../../$log
done
column -t ${m}-mer_info.txt >  "${m}-mer_info.txt.tmp" && mv "${m}-mer_info.txt.tmp" "${m}-mer_info.txt"
cd ..
cd ../../
footer wrote builder pdbs


#make other pdbs from monomeric pdbs according to other lists
header building other pdbs
cd pdb_bank/equal_atoms
for mer in ${compare_mer[@]}; do
    
    #dont do this part for 1-mers
    if [[ " ${mer} " == " 1 " ]]; then
        continue
    fi

    cd $mer-mer
    printf 'pdb-id\tatoms\tresidues\tchains\n' > ${mer}-mer_info.txt
    rm -f pdb_list_${mer}-mer.dat
    u_mer=$(wc -l ../../cut_and_cleaned/withoutgaps/${mer}-mer/pdb_list_${mer}-mer.dat | awk '{ print $1 }')
    for row in $(seq 1 "$u_mer") ; do
        id=$(awk "FNR == $row { print \$1 }" ../../cut_and_cleaned/withoutgaps/${mer}-mer/pdb_list_${mer}-mer.dat)
        
        #get chains to combine
        chainarr=()
        chains=${id: -${mer}}
        printf '\n'
        echo "${id}_withgaps" 2>&1 | tee -a ../../../$log
        while read -n1 char; do
            chainarr+=($char)
        done < <(echo -n "$chains")
        
        mkdir -p combining_here/
        newid="${id}_withgaps"
        newfile="${newid}.pdb"
        
        #copy to combine dir if they exist. if they dont, skip combo
        for c in ${chainarr[@]} ; do
            chainless_id="${id:0:$((${#id} - ${mer}))}" #removes last $mer characters (ie chain ids) from id    
            c_id="${chainless_id}${c}"
            c_file="${c_id}_withgaps.pdb"
            if [ ! -f "../1-mer/${c_file}" ]; then
                printf '\t' 2>&1 | tee -a ../../../$log
                echo "couldnt find ${c_file}. this is normal if you removed ${c_id} from the 1-mer list. moving on to next combo" 2>&1 | tee -a ../../../$log
                continue 2 
            fi
            cp ../1-mer/${c_file} combining_here/
            
            #combine one at a time
            cd combining_here
            cat $c_file >> $newfile
            printf '\t' 2>&1 | tee -a ../../../../$log
            echo "wrote $c_file into $newfile" 2>&1 | tee -a ../../../../$log
            cd .. 
        done
        mv combining_here/$newfile .
        rm -rf combining_here/
        printf '\t' 2>&1 | tee -a ../../../$log
        echo "reordering atoms and residues" 2>&1 | tee -a ../../../$log
        #renumber atoms
        python3 ../../../../code/renumber_atoms.py 1 ${newid}.pdb > ${newid}_a.pdb
        rm -f ${newid}.pdb
        #renumber residues
        python3 ../../../../code/renumber_res.py 1 ${newid}_a.pdb > ${newid}.pdb
        rm -f ${newid}_a.pdb
        

        #write info
        a=$(pdb_wc -a ${newid}.pdb | awk "FNR == 1 { print \$3 }")
        r=$(pdb_wc -r ${newid}.pdb | awk "FNR == 1 { print \$3 }")
        c=$(pdb_wc -c ${newid}.pdb | awk "FNR == 1 { print \$3 }")
        printf '%s\t%s\t%s\t%s\n' "${newid}" "$a" "$r" "$c" >> ${mer}-mer_info.txt

    done
    column -t ${mer}-mer_info.txt >  "${mer}-mer_info.txt.tmp" && mv "${mer}-mer_info.txt.tmp" "${mer}-mer_info.txt"
    cd ..
done
printf '\n' 2>&1 | tee -a $log

cd ../..

footer built pdbs
echo "the pdbs found in ${label}/pdb_bank/equal_atoms should be svd-ready, but there still might be unwanted pdbs." 2>&1 | tee -a $log
echo "files ${label}/pdb_bank/equal_atoms/*/*-mer_info.txt contain the number of atoms, residues, and chains for each pdb." 2>&1 | tee -a $log
echo "id check that to make sure its what you expect." 2>&1 | tee -a $log
printf '\n'

(printf '%0.s^v' {1..50}; printf '\n')  >> $log
(printf '%0.sv^' {1..50}; printf '\n')  >> $log
cd ..
