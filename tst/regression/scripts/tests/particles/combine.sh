#!/bin/bash

# Define the base name.
BASE=UniStream

# Find the number of blocks.
read nblocks < <(ls $BASE.block*.out1.00000.tab | wc -l)

# Find the number of snapshots.
read ntimes < <(ls $BASE.block0.out1.*.tab | wc -l)
(( --ntimes ))

for i in $(seq -f "%05g" 0 $((ntimes - 1))); do
	dst=$BASE.pout.$i.tab
	if (( nblocks > 1 )); then
		grep --regexp="^ *#" $BASE.block0.pout.$i.tab > $dst
		grep --regexp="^ *#" --invert-match --no-filename \
				$BASE.block*.pout.$i.tab | \
			sort --general-numeric-sort >> $dst
	else
		mv $BASE.block0.pout.$i.tab $dst
	fi
	cat $BASE.block*.out1.$i.tab > $BASE.out1.$i.tab
done

find . -name $BASE'.block*.*.tab' -delete
