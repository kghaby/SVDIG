#!/bin/bash
set -ae 

#STEPS
#==================================================
#step1.sh
	#blast, output plot of e values
#step2.sh
	#download pdbs, align, output alignment of pdb sequences
#step3.sh
	#cuts pdbs and separates into mers
#==================================================

#VARIABLES
#==================================================

#RELEVANT TO STEP 1 (and subsequent steps)
#=========================
label="mcm_test"
        #eg. label="mcm" will make dir mcm
        #use directory-appropriate formatting. also dont use spaces

fasta="fasta_bank/rcsb_pdb_3JA8_mcm2.fasta"
        #eg. fasta_bank/rcsb_pdb_3JA8_mcm2.fasta will be blasted
        #can only handle one sequence
        #need to be a file
ethresh="10"
        #eg. ethresh="10" will gather pdb's with e values under 10
        #can use python scientific notation form too, eg "1e-3" instead of "0.001"
notes="
                - https://www.youtube.com/watch?v=dQw4w9WgXcQ"
        #eg. notes=" - used 3ja8 sequence" will append this note to the log corresponding to $label
#=========================

#RELEVANT TO STEP 2 (and subsequent steps)
#=========================
newethresh="1e-60"
        #eg. ethresh="10" will use pdb's with e values under 10
        #can use python scientific notation form too, eg "1e-3" instead of "0.001"
pdb_ref_id="3JA8"
#use the id as it is named in the <label>/pdb_bank/alignment directory (minus the ".pdb") in case the original pdb had multiple models, lettered and numbered chain id's, etc
#pdb_ref_chains=(A B C D E F)
    #these are the chains in ref that you want to standardize across all pdbs
    #spaces between chains
    #script uses all chains
#=========================

#RELEVANT TO STEP 3
#=========================
residue_range=(900:1060 1720:1855)
        #positions of residues in the alignment output (eg alignment_output.fasta)
        #eg residue_range=(900:1060 1720:1855) 
compare_mer=(1 2)
        #the oligomerization state(s) you want to compare in svd
        #this script supports these types of oligomers: compare_mer=(1 2 4 6)
                #easy to add more types if you want.
                        #1.  add an nmer_pattern array with the other nmer_pattern arrays
                        #2.  down near the bottom of the script, under the other "splitpdb ..." lines, add another splitpdb line for the n-mer of your choice. change the 1st and 2nd arguments, but keep the format
        #eg compare_mer=(2) for dimer-dimer comparisons, compare_mer=(1 2) for monomer-monomer and dimer-dimer comparisons
        #indep of seq ie can use a monomer seq but compare dimers
dimer_pattern=(A,D D,B B,F F,C C,E E,A)
tetramer_pattern=(A,B,C,D E,F,G,H)
hexamer_pattern=(A,B,C,D,E,F)
        #only need to edit the >1-mers listed in $compare_mer
        #must be capital letters
        #this is basically a set of the chain combos you want to compare in svd
        #if the sequence aligns to a pdb with less chain thans the reference, this script will make whatever combos possible out of the elements in these oligo patterns
        #eg dimer_pattern=(A,D D,B B,F F,C C,E E,A) creates pdbs for each interacting pair available in a hexameric MCM where the chains are ordered around the circle as ADBFCE
        #eg hexamer_pattern=(A,B,C,D,E,F) creates pdbs containing chains A-F of downloaded pdb's
#=========================

#TBD
#=========================
#atoms="CA,C,N,O"
        #eg. atoms="CA,C,N,O" will select backbone atoms from pdb
        #chosen later

#alltitles="no"
        #eg. alltitles="no" will only extract one pdb per matching sequence
        #has not been implemented yet
#alignment="0"
        #"0" to use the alignment from blast
        #"1" to run another alignment
#==================================================

#FUNCTIONS
#==================================================
header () {
        printf '\n%s\n' "$*..." 2>&1 | tee -a $log
        (printf '%0.s=' {1..80}; printf '\n') 2>&1 | tee -a $log
}

footer () {
        (printf '%0.s=' {1..80}; printf '\n') 2>&1 | tee -a $log
        printf '%s\n\n' "...$*" 2>&1 | tee -a $log
}

checkchains () {
        checkchn=$(pdb_wc $1 | awk "FNR == 2 { print \$3 }")
        if [ ! $checkchn -eq $2 ] ; then
                rm $1
                continue
        fi
}

splitpdb () {
        if [[ " ${compare_mer[@]} " =~ " $1 " ]]; then
                local mer="$1" # Save first argument in a variable
                shift   # Shift all arguments to the left (original $1 gets lost)
                pattern=("$@") # Rebuild the array with rest of arguments
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
                for chn in $(seq 0 "$((${#pattern[@]}-1))"); do
                        set +e
                        nocomma=$(echo ${pattern[$chn]} | tr -d ',')
                        pdb_selchain -${pattern[$chn]} $id.pdb > ${id}_$nocomma.pdb
                        checkchains ${id}_$nocomma.pdb $mer
                        echo "    ${id}_$nocomma.pdb" 2>&1 | tee -a ../../../$log
                        printf "yes,${id}_$nocomma,$nocomma,,$nocomma \n" >> ../../../metadata/${label}_svd_$mer-mer.csv
                        set -e
                done
                rm $id.pdb
                cd ..
        fi
}

skip_if_not_found () {
        if [ ! -f "$1" ]; then
                echo "$1 was not found. moving on to the next one" 2>&1 | tee -a ../../$log
                continue
        fi
}

#usage: print_array "${array_name[@]}"
print_array () 
{
  # run through array and print each entry:
  local array
  array=("$@")
  for i in "${array[@]}" ; do
      printf '%s ' "$i"
  done
  # print a new line
  printf '\n'
}

#==================================================
