#!/bin/bash

datasets=("yh3")
#datasets=("yh1")
nreps=100
path="."
qspec="-q normal"

for ds in "${datasets[@]}"; do
	for ((i=1;$i<=$nreps;i++)); do
		prefix="$path/$ds/data.$i.orig.boot"
		cmd="./mpboot-avx -s $path/$ds/data.$i -bb 1000 -keep_ident -pure_sa -cooling_schedule 0 -pre $prefix -seed 213 > $prefix.out 2>&1"
		eval $cmd
        # bsub -J bo$ds$i $qspec $cmd
	done
done


