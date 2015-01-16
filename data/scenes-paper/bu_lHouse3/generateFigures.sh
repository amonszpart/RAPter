#!/bin/sh

SCRIPT_PATH="../../../globOpt/scripts"

mkdir figure

############################################################
## generate stats from input

#distribution
python ../normal_distr.py cloud.ply ndistrCloudBinary.svg "House - Input"

meshlabserver -i cloud.ply -o figure/cloud_pts.ply -s $SCRIPT_PATH/colorize.mlx -om vc vn

splatting figure/cloud_pts.ply figure/cloud.ply 2 0.015


############################################################
## generate stats for patches
../globOptVis  --show3D  --pop-limit 3 -p patches.csv -a points_primitives.csv --title GlobOpt - [Dir-Colours] patches --angle-gens 0 --draw-mode 28 --save-poly  --paral-colours 

#distribution
python ../normal_distr.py cloudRGBNormal_patches_reProj_noUnass_noPrim.ply ndistrPatches.svg "House - Patches"


############################################################
## generate cleaned point-cloud model, and associated stats
../globOptVis  --show3D  --pop-limit 3 -p primitives_it20.bonmin.csv -a points_primitives_it19.csv --title GlobOpt - [Dir-Colours] 20 iteration output --angle-gens 0 --draw-mode 28 --save-poly  --paral-colours 

#distribution
python ../normal_distr.py cloudRGBNormal_it20_reProj_noUnass_noPrim.ply ndistrIt20.svg "House - Iteration 20"

meshlabserver -i cloudRGBNormal_it20_reProj_noUnass_noPrim.ply -o figure/cloud_cleaned_pts.ply -om vc vn

splatting figure/cloud_cleaned_pts.ply figure/cloud_cleaned.ply 2 0.015

#graph
python ../readGraphProperties.py primitives_it20.bonmin.csv points_primitives_it19.csv cloud.ply --angles 0 --iteration 20

############################################################
## generate planar approximation
../globOptVis  --show3D  --pop-limit 3 -p primitives_it20.bonmin.csv -a points_primitives_it19.csv --title GlobOpt - [Dir-Colours] 10 iteration output --angle-gens 0 --draw-mode 1 --save-poly  --paral-colours --no-pts
