#!/bin/bash

echo "--TSP Benchmark--"

make

for instance in instances/*; do
	echo $instance >> ./benchmark/bm.txt
	for i in {1..10}; do
		./tsp ${instance} | grep 'COST\|TIME' | awk "{print $1}" >> ./benchmark/bm.txt
	done

done

echo "-" >> ./benchmark/bm.txt

cd benchmark
python3 bm.py
cd ..

echo "Benchmark completed."