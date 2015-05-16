#!/bin/bash

function print_usage() {
        echo "usage:\t showcands.sh iteration [3D] [dids]" 
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
	if [[ ! $2 == "3D" ]]; then
		dids=$2
	else
		flag3D="3D"
		if [[ ! -z "$3" ]]; then
			dids="--dids $3"
		fi
	fi
fi

prims="candidates_it$iteration.csv";
assignments="points_primitives_it$prev.csv";
echo "$prims $assignments"
../globOptVis --show$flag3D --scale 0.1 --pop-limit 2 --title "GlobOpt - $iteration iteration candidates ($dids)" --angle-gens 90 --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-scale --bg-colour .9,.9,.9 --no-rel $dids -p $prims -a $assignments &

echo "../globOptVis --show$flag3D --scale 0.1 --pop-limit 2 --title \"GlobOpt - $iteration iteration candidates ($dids)\" --angle-gens 90 --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-scale --bg-colour .9,.9,.9 --no-rel $dids -p $prims -a $assignments &"



# cone + missing
# angular plots: 1. perfect solution, 2. initial lines, 3. solution (2.: manual correpsondance)
# thetaphi diagram : initial many, in final only the correct number (don't need to colour for now)
# paper: initial results section -> tell Niloy
