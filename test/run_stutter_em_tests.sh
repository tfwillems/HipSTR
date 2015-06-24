count=1
for nsamp in 100 250 500
do
    for mu in 0.01 0.001 0.0001
    do
	for beta in 0 0.25 0.5
	do
	    for p_geom in 0.75 0.9 1.0
	    do
		for stutter_geom in 0.8 0.95
		do
		    for stutter_up in 0.01 0.15
		    do
			for stutter_down in 0.01 0.15
			do
			    for read_counts in 1,2,3 2,3,4 3,4,5 10,11,12
			    do
				for iter in $(seq 1 25)
				do
				    echo $nsamp $mu $beta $p_geom $stutter_geom $stutter_up $stutter_down $read_counts $count
				    count=`expr $count + 1`
				done
			    done
			done
		    done
		done
	    done
	done
    done
done | xargs -L 1 -P 20 ./run_stutter_em_test.sh
