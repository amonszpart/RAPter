#!/bin/bash
function printUsage()
{
	echo "Usage: \"./showSchnabel name [3D]\" to display primitives.$1.csv"
}

if [[ -z "$1" ]]; then printUsage; exit 1; else name=$1; fi
if [[ -z "$2" ]]; then flag3D=""; else flag3D=$2; fi

echo "../globOptVis --show$flag3D --scale 0.01 -p $name.primitives.csv -a $name.points_primitives.csv --cloud $name.cloud.ply --pop-limit 3 --title \"Ransac $name output\" --angle-gens 90 --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-scale --bg-colour 1,1,1 --no-rel --paral-colours --perfect-angle 0.0001 &"
../globOptVis --show$flag3D --scale 0.01 -p $name.primitives.csv -a $name.points_primitives.csv --cloud $name.cloud.ply --pop-limit 3 --title \"Ransac $name output\" --angle-gens 90 --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-scale --bg-colour 1,1,1 --no-rel --paral-colours --perfect-angle 0.0001 &