#!/bin/bash
if [[ -z "$1" ]]; then it=20; else it=$1; fi

rm -r *.bak representatives*.csv *.gv candidates_representatives* *lvl1* problem points_representatives*.csv *rprtmp*
itm1=$(( $it - 1 ))
echo "git add segments.csv points_segments.csv patches.csv points_primitives.csv primitives_it${it}.bonmin.csv points_primitives_it$itm1.csv run.log"