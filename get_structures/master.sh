#!/bin/bash
set -ae 

#EXAMPLE WORKFLOW
#==================================================
#step1.sh
    	#blast, output plot of e values
#step2.sh
    	#download pdbs, separate models and numbered/lettered chains
#step3.sh
    	#standardize chains, get sequences, align sequences
#step4.sh
    	#cut pdbs and separate into n-mers
#step5.sh 
    	#fill gaps to make all pdbs have the same number of residues (and atoms if only backbone atoms are kept)
    	#only keeps atoms specified
#==================================================

#VARIABLES
#==================================================

#RELEVANT TO ALL STEPS
#=========================
label="example"
        #eg. label="example" will make dir example
        #use directory-appropriate formatting. also dont use spaces
#=========================

#RELEVANT TO STEP 1
#=========================
fasta="fasta_bank/example.fasta"
        #eg. fasta="fasta_bank/example.fasta"
        #can only handle one sequence per fasta
        #needs to be a fasta file
ethresh="10"
        #eg. ethresh="10" will gather pdb's with e values under 10
        #can use python scientific notation form too, eg "1e-3" instead of "0.001"
notes="
                - https://www.youtube.com/watch?v=dQw4w9WgXcQ"
        #eg. notes=" - used 3ja8 sequence" will append this note to the log corresponding to $label
#=========================

#RELEVANT TO STEP 2
#=========================
newethresh="1e-3"
        #eg. ethresh="1" 
            #will use pdb's with e values under 1
        #can use python scientific notation form too, eg "1e-3" instead of "0.001"
#=========================

#RELEVANT TO STEP 3
#=========================
pdb_ref_id="3JA8"
    #this is the name of the pdb file you wish to use, ie the name of the file minus ".pdb". 
    #assuming youve done step 2, check $label/pdb_bank/aligned/pdb_list_standard_input.dat for the list of available pdb's. note that the name would be different from the four digit pdb id if there contained multiple models or lettered and numbered chain id's
        #this is also an opportunity to manually remove from pdbs by either deleting them from pdb_bank/aligned or deleting the pdb's row in pdb_list_standard_input.dat  
    #eg. pdb_ref_id="6M0J_ltr"
        #will use 6M0J_ltr.pdb
#=========================

#RELEVANT TO STEP 4
#=========================
alignment_file="${label}/pdb_bank/aligned/alignment_files/align_output.fasta"
    #this is the output of step 3 alignment that will be used to fill gaps in the full structure
    #if you did your own alignment, this should point to the output fasta file with a path relative to the main directory with the step*.sh files
    #eg alignment_file="${label}/pdb_bank/aligned/alignment_files/align_output.fasta"
        #uses the alignment output from step 3
residue_range=(904:1100 1715:1956 2608:2703 3495:3633 4189:4289 4891:5073)
        #positions of residues in the alignment output (eg alignment_output.fasta)
        #eg residue_range=(900:1060 1720:1855) 
compare_mer=(1 2 6)
        #the oligomerization state(s) you want to compare in svd
        #this script supports these types of oligomers: compare_mer=(1 2 4 6)
                #easy to add more types if you want.
                        #1.  add an nmer_pattern array with the other nmer_pattern arrays
                        #2.  down near the bottom of the script, under the other "splitpdb ..." lines, add another splitpdb line for the n-mer of your choice. change the 1st and 2nd arguments, but keep the format
        #eg. compare_mer=(2) 
            #for dimer-dimer comparisons, 
        #eg. compare_mer=(1 2) 
            #for monomer-monomer and dimer-dimer comparisons
        #indep of seq ie can use a monomer seq but compare dimers
        #1 is done no matter what since the monomers are used to build the larger pdbs, so you dont need to include 1.
        #also used in step 5
dimer_pattern=(A,D D,B B,F F,C C,E E,A)
tetramer_pattern=(A,B,C,D E,F,G,H)
hexamer_pattern=(A,B,C,D,E,F)
        #assume first chain is A 
        #only need to edit the >1-mers listed in $compare_mer
        #must be capital letters
        #this is basically a set of the chain combos you want to compare in svd
        #if the sequence aligns to a pdb with less chain thans the reference, this script will make whatever combos possible out of the elements in these oligo patterns
        #eg. dimer_pattern=(A,D D,B B,F F,C C,E E,A) 
            #creates pdbs for each interacting pair available in a hexameric MCM where the chains are ordered around the circle as ADBFCE
        #eg. hexamer_pattern=(A,B,C,D,E,F) 
            #creates pdbs containing chains A-F of downloaded pdb's
#=========================

#RELEVANT TO STEP 5
#=========================
compare_mer=$compare_mer
    #compare_mer from step 4
atoms=(CA C N O)
    #eg. atoms=(CA C N O) 
        #will select backbone atoms from pdb
#IMPORTANT NOTE
    #the list of pdbs that will combine to form the rest of the oligo states is "${label}/pdb_bank/cut_and_cleaned/withoutgaps/1-mer/pdb_list_1-mer.dat" 
        #this is a list of monomeric pdbs
        #edit this list and delete any lines of the unwanted pdbs 
        #this is recommended as some unwanted pdbs may be in there, like random chains of dna. this may mess up the alignment, and will add unnecessary gaps
#=========================

#
#
#
#
#
#
#
#
#
#
#
#
#
#
#

#TO DO
#=========================
#alltitles="no"
        #eg. alltitles="no" will only extract one pdb per matching sequence
        #has not been implemented yet
#pdb_ref_chains=(A B C D E F)
    #these are the chains in ref that you want to standardize across all pdbs
    #spaces between chains
    #script uses all chains
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
                rm -f $1
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
                echo "separating $id into $mer-mers..." 2>&1 | tee -a ../../../../$log
                if [ $chains -lt $mer ] ; then
                        echo "$id has less than $mer chains (it has $chains) so it cant be used for $mer-mer comparisons." 2>&1 | tee -a ../../../../$log
                        rm -f $id.pdb
                        cd ..
                        continue
                fi
                echo "from chains $allchains, wrote:" 2>&1 | tee -a ../../../../$log
                for chn in $(seq 0 "$((${#pattern[@]}-1))"); do
                        set +e
                        nocomma=$(echo ${pattern[$chn]} | tr -d ',')
                        python3 ../../../../../code/take_chains.py -${pattern[$chn]} $id.pdb > ${id}_$nocomma.pdb

                        #if the last line is not TER, add it
                        F="${id}_${nocomma}.pdb"
                        if [[ ! $(tail -n 1 $F | grep 'TER') ]]; then
                            (printf "TER"; printf '%0.s ' {1..77}; printf '\n') >> $F
                        fi
 
                        #if the first line is TER, remove it
                        F="${id}_${nocomma}.pdb"
                        if [[ $(head -n 1 $F | grep 'TER') ]]; then
                            tail -n +2 "$F" >  "$F.tmp" && mv "$F.tmp" "$F"
                        fi

                        checkchains ${id}_$nocomma.pdb $mer
                        echo "${id}_${nocomma}" >> pdb_list_${mer}-mer.dat
                        echo "    ${id}_$nocomma.pdb" 2>&1 | tee -a ../../../../$log
                        printf "yes,${id}_$nocomma,$nocomma,,$nocomma \n" >> ../../../../metadata/${label}_svd_$mer-mer.csv
                        set -e
                done
                rm -f $id.pdb
                cd ..
        fi
}

skip_if_not_found () {
        if [ ! -f "$1" ]; then
                echo "$1 was not found. moving on to the next one" 2>&1 | tee -a ../../$log
                continue
        fi
}


ter_first () {
    #if the first line is TER, remove it
    F=$1
    if [[ $(head -n 1 $F | grep 'TER') ]]; then #$(sed -n '/^TER/p;q' $F) is faster but maybe not portable?
        tail -n +2 "$F" >  "$F.tmp" && mv "$F.tmp" "$F"
    fi
}

ter_last () {
    #if the last line is not TER, add it
    F=$1
    if [[ ! $(tail -n 1 $F | grep 'TER') ]]; then
        (printf "TER"; printf '%0.s ' {1..77}; printf '\n') >> $F
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
