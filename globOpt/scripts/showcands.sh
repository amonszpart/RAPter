#!/bin/bash

function print_usage() {
        echo "usage:\t showcands.sh iteration" 
}

# parse iteration
if [[ -z "$1" ]]; then
	print_usage
	exit 1
else
	iteration=$1;
	prev=$(($iteration-1))
fi

if [[ ! -z "$2" ]]; then
	dids=$2
fi

prims="candidates_it$iteration.csv";
assignments="points_primitives_it$prev.csv";
echo "$prims $assignments"
../globOptVis --show --scale 0.019 --pop-limit 5 -p $prims -a $assignments --title "GlobOpt - $iteration iteration candidates ($dids)" --angle-gens 60,90 --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-scale --bg-colour .9,.9,.9 --ids --no-rel --dids $dids &

echo "../globOptVis --show --scale 0.019 --pop-limit 5 -p $prims -a $assignments --title \"GlobOpt - $iteration iteration candidates\" --angle-gens 60,90 --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-scale --bg-colour .9,.9,.9 --ids --no-rel"



# cone + missing
# angular plots: 1. perfect solution, 2. initial lines, 3. solution (2.: manual correpsondance)
# thetaphi diagram : initial many, in final only the correct number (don't need to colour for now)
# paper: initial results section -> tell Niloy
