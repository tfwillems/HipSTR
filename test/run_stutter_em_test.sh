#!/bin/bash

nsamp=$1
mu=$2
beta=$3
p_geom=$4
stutter_geom=$5
stutter_up=$6
stutter_down=$7
read_counts=$8
count=$9

# Simulate the STRs and the set of observed reads
python ../../str_mut_models/main.py --mu $mu --beta $beta --pgeom $p_geom --sim_str_reads --tree ../../str_mut_models/data/1000Y.MLtree.lengths.bootstraps.rooted.rotated.20140826.nwk --out garb$count --nsamp $nsamp --stutter_inc $stutter_up --stutter_dec $stutter_down --stutter_geom $stutter_geom --obs_read_counts $read_counts

for haploid in False True
do
    if [ $haploid = "True" ]
    then
	freqs=(0.01)
    else
	freqs=(0.01 0.25 0.75)
    fi


    for phase_freq in ${freqs[@]}
    do	
	# Estimate the underlying stutter model
	res=`./em_stutter_test garb$count.str_reads.txt $haploid 4 $phase_freq | cut -f 4-`
	echo $mu $beta $p_geom $stutter_geom $stutter_down $stutter_up $haploid $nsamp $read_counts $phase_freq  $res >> stutter_results.txt
    done
done

# Delete the file of read counts
rm garb$count.str_reads.txt
