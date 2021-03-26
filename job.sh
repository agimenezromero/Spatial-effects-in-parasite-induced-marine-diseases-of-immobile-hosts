#!/bin/bash

export LC_NUMERIC="en_US.UTF-8"

for param in $(seq 0.0 0.02 10) 
do	
	#run -c 10 -o kappa_$param.out -e kappa_$param.err julia Cell_IBM_Nuredduna.jl "$param"
    run -c 10 -m 70 julia Cell_IBM_Nuredduna.jl "$param"
	#run -x -c 10 -m 70 julia Cell_IBM_Nuredduna.jl "$param"
	
	sleep 0.05
done
