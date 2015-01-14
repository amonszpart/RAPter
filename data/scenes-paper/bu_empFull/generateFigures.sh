#!/bin/sh

SCRIPT_PATH="../../../globOpt/scripts"

ln -s cloudBinary.ply cloud.ply
mkdir figure

############################################################
## generate stats from input

#distribution
python ../normal_distr.py cloud.ply ndistrCloudBinary.svg "Euler - Input"

meshlabserver -i cloud.ply -o figure/cloud.ply -s $SCRIPT_PATH/colorize.mlx -om vc vn

splatting figure/cloud.ply figure/cloud.ply 0 0.001


############################################################
## generate stats for patches
../globOptVis  --show3D  --pop-limit 3 -p patches.csv -a points_primitives.csv --title GlobOpt - [Dir-Colours] patches --angle-gens 0 --draw-mode 28 --save-poly  --paral-colours 

#distribution
python ../normal_distr.py cloudRGBNormal_patches_reProj_noUnass_noPrim.ply ndistrPatches.svg "Euler - Patches"

meshlabserver -i cloudRGBNormal_patches_reProj_noUnass_noPrim.ply -o figure/cloud_cleaned.ply -om vc vn

splatting figure/cloud_cleaned.ply figure/cloud_cleaned.ply 0 0.001


############################################################
## generate cleaned point-cloud model, and associated stats
../globOptVis  --show3D  --pop-limit 3 -p primitives_it10.bonmin.csv -a points_primitives_it9.csv --title GlobOpt - [Dir-Colours] 10 iteration output --angle-gens 0 --draw-mode 28 --save-poly  --paral-colours 

#distribution
python ../normal_distr.py cloudRGBNormal_it10_reProj_noUnass_noPrim.ply ndistrIt69.svg "Euler - Iteration 69"

meshlabserver -i cloudRGBNormal_it10_reProj_noUnass_noPrim.ply -o figure/cloud_cleaned_pts.ply -om vc vn

splatting figure/cloud_cleaned_pts.ply figure/cloud_cleaned.ply 0 0.002

#graph
python ../readGraphProperties.py primitives_it10.bonmin.csv points_primitives_it09.csv cloud.ply --angles 0 --iteration 10

############################################################
## generate planar approximation
../globOptVis  --show3D  --pop-limit 3 -p primitives_it10.bonmin.csv -a points_primitives_it9.csv --title GlobOpt - [Dir-Colours] 10 iteration output --angle-gens 0 --draw-mode 28 --save-poly  --paral-colours --no-pts
