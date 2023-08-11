#!/bin/bash


nreps=100
datasets=("yh3")
#datasets=($2)
path="."

for ds in "${datasets[@]}"; do
	for ((i=1;$i<=$nreps;i++)); do
		prefix="data.$i.orig.boot"
		cmd="./mpboot-avx -rf tree.$i $prefix.treefile -v"
                eval $cmd
	done

done



