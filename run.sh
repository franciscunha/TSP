#!/bin/bash

echo "Running TSP"

make

for instance in ./instances/*; do
	echo $instance
	for i in {1..10}; do
		./tsp ${instance} | grep 'COST\|TIME' | awk "{print $1}"
	done

done
