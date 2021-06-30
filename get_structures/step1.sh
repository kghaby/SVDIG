#!/bin/bash
set -ea; trap 'echo "Mission failed, we'"'"'ll get line $LINENO next time."' ERR
source master.sh

#MAIN
#==================================================

#check if label has been used
if [ -d "$label" ]; then
	echo "The label $label has already been used. Continuing will overwrite: 
	blast_results.xml
	any pdbs that have the same name as those found in this run
	alignment data files
	evalues_$label.dat
	evalues_$label.png"
	if [ -f "$label/blast_files/hits_ethresh$ethresh.dat" ]; then
		echo "	hits_ethresh$ethresh.dat"
	fi
	csv_ar=(mcm/metadata/*csv)
	for csv in ${csv_ar[@]}; do
		if [ -f "$csv" ]; then
			echo "	$csv"
		fi	
	done
	select yn in "Use existing directory" "Exit script"; do
		case $yn in
			"Use existing directory" ) echo "Moving on..."; break;;
			"Exit script" ) echo "Change label at the top of this file"; exit;;
		esac
	done
fi


#make directories
mkdir -p $label/{pdb_bank/{crude,aligned,cut_and_cleaned},blast_files,metadata}
cd $label
fasta="../$fasta"

#write log
log="$label".log
printf '\n%s' "Step 1 for $label - " 2>&1 | tee -a $log
date 2>&1 | tee -a $log
printf 'Step 1 Variables\n\tNotes: %s\n\tE value thresh for blast: %s\n\tAtoms: %s\n\tFasta: %s\n\tFasta Contents: %s\n\n' \
	"$notes" "$ethresh" "$atoms" "$fasta" "$(<$fasta)" 2>&1 | tee -a $log

#append a line to README for convenience
printf "$label - " >> README
date >> README


#search blast
header blasting
echo "will take a sec"
python3 ../code/blast.py 2>&1 | tee -a $log
echo "Done with code/blast.py" 2>&1 | tee -a $log
python3 ../code/compare_e.py > blast_files/evalues_$label.dat
echo "Done with code/compare_e.py. wrote blast_files/evalues_$label.dat" 2>&1 | tee -a $log
python3 ../code/parse_blast.py > blast_files/hits_ethresh$ethresh.dat
echo "Done with code/parse_blast.py. wrote blast_files/hits_ethresh$ethresh.dat" 2>&1 | tee -a $log
cp $fasta blast_files
footer blasted


#generate threshold plot 
header plotting
num_eval=$(($(wc -l blast_files/evalues_$label.dat | awk '{ print $1 }')-1))
echo "amount of e values in $label/blast_files/evalues_$label.dat: $num_eval" 2>&1 | tee -a $log 
fontsize=$((500 / $num_eval))
xres=1000
if [ $fontsize -ge 11 ]; then
	echo "changing fontsize from $fontsize to 10" 2>&1 | tee -a $log
	fontsize=10
fi
if [ $fontsize -le 7 ]; then
	echo "changing fontsize from $fontsize to 8." 2>&1 | tee -a $log
	newxres=$(($((8-$fontsize)) * 200 + $xres))
	fontsize=8
	echo "	To compensate, changing plot length from $xres to $newxres" 2>&1 | tee -a $log
	xres=$newxres
fi
gnuplot -persist << EOF
	set terminal png size $xres,600
	set output 'blast_files/evalues_$label.png'
	set multiplot layout 2,1 rowsfirst
	TMARGIN = "set tmargin at screen 0.90; set bmargin at screen 0.54"
	BMARGIN = "set tmargin at screen 0.54; set bmargin at screen 0.18"
	SIDEMARGIN = "set lmargin at screen 0.15; set rmargin at screen 0.95"
	set title "Unique E values of $label BLAST hits"
	set ytics offset 0,0.5
	unset xtics
	@TMARGIN; @SIDEMARGIN
	plot "blast_files/evalues_$label.dat" u 2:xtic(1) w lp lt 1 pt 31 notitle
	set notitle
	set logscale y 10
	set ytics offset 0,-0.5
	set xtics rotate by -45 font ", $fontsize"
	set xlabel "PDB|CHAIN (1 per e value)"
	set label 1 "E value" at screen 0.03,0.54 center rotate by 90
	@BMARGIN; @SIDEMARGIN
	plot "blast_files/evalues_$label.dat" u 2:xtic(1) w lp lt 1 pt 31 notitle
EOF
set +e;feh -qF blast_files/evalues_$label.png & echo "$label/blast_files/evalues_$label.png should be displayed if you have feh" 
set -e
footer plotted


(printf '%0.s^v' {1..50}; printf '\n')  >> $log
(printf '%0.sv^' {1..50}; printf '\n')  >> $log
cd ..
