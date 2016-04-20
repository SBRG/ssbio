#!/bin/bash

for i in *.pdb; do
    echo "Running ${i}..."
    seq=`echo ${i}  | cut -d. -f1`
    curr_dir=`pwd`

	cp ${i} temp.pdb

    tleap -f leaprc

    rm temp.pdb
	mv temp_mod.pdb xleap_modified/${seq}_xleap.pdb
	mv temp_mod.prmtop amber_minimized/${seq}.prmtop
	mv temp_mod.inpcrd amber_minimized/${seq}.inpcrd

done
