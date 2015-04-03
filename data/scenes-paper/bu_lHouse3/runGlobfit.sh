#!/bin/bash
scale=0.005
poplimit=3
anglegens="90"
visdefparam="--angle-gens $anglegens --use-tags --no-clusters --statuses -1,1 --no-pop --no-rel --no-scale --bg-colour .9,.9,.9" #"--use-tags --no-clusters" #--ids

export PATH="/home/bontius/matlab_symlinks/":$PATH
export WORKSPACE=/export/home/kandinsky/nmellado/workspace_globOpt
$WORKSPACE/globOpt/globOpt/build/Release/bin/toGlobFit --planes --prims segments.csv --cloud cloud.ply -a points_segments.csv --scale $scale
$WORKSPACE/3rdparty/globfit/build/bin/globfit_release -i segments.globfit -v -o 3.0 -g 3.0 -a 0.1 -p $scale -l $scale -r $scale
$WORKSPACE/globOpt/globOpt/build/Release/bin/toGlobFit --from segments_oa.globfit --planes --prims segments.csv --cloud cloud.ply -a points_segments.csv --scale $scale

../../globOptVis --show3D --scale $scale --pop-limit $poplimit -p primitives_it9.bonmin.csv -a points_primitives_it8.csv --title "GlobOpt" $visdefparam --paral-colours --no-rel &

../../globOptVis --show3D --scale $scale --pop-limit $poplimit -p primitives.globfit.csv -a points_primitives.globfit.csv --title "Globfit" $visdefparam --paral-colours --no-rel &

../../globOptVis --show3D --scale $scale --pop-limit $poplimit -p patches.csv -a points_primitives.csv --title "Input" $visdefparam --paral-colours --no-rel &
