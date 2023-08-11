#!/bin/bash

dir="."
dslist=(
	# 'yh1'
	# 'yh2'
	'yh3'
	# 'yh4'
	# 'yh5'
	)

type=(
	# 'dna'
	# 'dna'
	'dna'
	# 'prot'
	# 'prot'
	)

ntaxa=(
	# '100'
	# '200'
	'500'
	# '100'
	# '200'
	)

nsites=(
	# '500'
	# '1000'
	'1000'
	# '300'
	# '500'
	)
		
nreps=100

count=0
while [ "x${dslist[count]}" != "x" ]
do	
	dataset=${dslist[count]}

	cd $dir
	mkdir $dataset	
	cd $dir/$dataset
	
	if [ "${type[count]}" == "dna" ]; then
		for ((i=1;i<=$nreps;i++)); do 
			../mpboot-avx -r ${ntaxa[count]} -rlen 0.02 0.1 0.3 -rzero 0 tree.$i
			../seq-gen -z$RANDOM -mHKY -t1.0 -l${nsites[count]} tree.$i > data.$i
			../mpboot-avx -s data.$i -af fasta -ao data.$i.fa
		done
	fi

	if [ "${type[count]}" == "prot" ]; then
		for ((i=1;i<=$nreps;i++)); do 
			../mpboot-avx -r ${ntaxa[count]} -rlen 0.02 0.1 0.3 -rzero 0 tree.$i
			../seq-gen -z$RANDOM -mWAG -l${nsites[count]} tree.$i > data.$i
			../mpboot-avx -s data.$i -af fasta -ao data.$i.fa
		done
	fi

	count=$(( $count + 1 ))
done

