#!/bin/bash

function print_usage() {
        echo "usage:\t showquick.sh iteration [3D]" 
}

# parse iteration
if [[ -z "$1" ]]; then
	print_usage
	exit 1
else
	iteration=$1;
	prev=$(($iteration-1))
fi

if [[ -z "$2" ]]; then flag3D=""; else flag3D="$2"; fi

prims="primitives_it$iteration.bonmin.csv";
assignments="points_primitives_it$prev.csv";
echo "$prims $assignments"
../globOptVis --show$flag3D --scale 0.02 --pop-limit 3 -p $prims -a $assignments --title "GlobOpt - $iteration iteration output" --angle-gens 60,90 --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-scale --bg-colour .9,.9,.9 --no-rel &

echo "../globOptVis --show$flag3D --scale 0.02 --pop-limit 3 --title \"GlobOpt - $iteration iteration output\" --angle-gens 90 --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-scale --bg-colour .9,.9,.9 --no-rel -p $prims -a $assignments"



# cone + missing
# angular plots: 1. perfect solution, 2. initial lines, 3. solution (2.: manual correpsondance)
# thetaphi diagram : initial many, in final only the correct number (don't need to colour for now)
# paper: initial results section -> tell Niloy
