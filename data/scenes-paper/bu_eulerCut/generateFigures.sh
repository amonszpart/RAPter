#!/bin/sh

SCRIPT_PATH="../../../globOpt/scripts"

mkdir figure

############################################################
## generate stats from input

#distribution
python ../normal_distr.py cloud.ply ndistrCloudBinary.svg "Euler - Input"

meshlabserver -i cloud.ply -o figure/cloud_pts.ply -s $SCRIPT_PATH/colorize.mlx -om vc vn

splatting figure/cloud.ply figure/cloud_pts.ply 1 0.001


############################################################
## generate stats for patches
../globOptVis  --show3D  --pop-limit 3 -p patches.csv -a points_primitives.csv --title GlobOpt - [Dir-Colours] patches --angle-gens 0 --draw-mode 28 --save-poly  --paral-colours 

#distribution
python ../normal_distr.py cloudRGBNormal_patches_reProj_noUnass_noPrim.ply ndistrPatches.svg "Euler - Patches"

meshlabserver -i cloudRGBNormal_patches_reProj_noUnass_noPrim.ply -o figure/cloud_cleaned.ply -om vc vn



############################################################
## generate cleaned point-cloud model, and associated stats
../globOptVis  --show3D  --pop-limit 3 -p primitives_it45.bonmin.csv -a points_primitives_it44.csv --title GlobOpt - [Dir-Colours] 45 iteration output --angle-gens 0 --draw-mode 28 --save-poly  --paral-colours 

#distribution
python ../normal_distr.py cloudRGBNormal_it45_reProj_noUnass_noPrim.ply ndistrIt45.svg "Euler - Iteration 45"

meshlabserver -i cloudRGBNormal_it45_reProj_noUnass_noPrim.ply -o figure/cloud_cleaned_pts.ply -om vc vn

splatting figure/cloud_cleaned_pts.ply figure/cloud_cleaned.ply 1 0.001

#graph
python ../readGraphProperties.py primitives_it45.bonmin.csv points_primitives_it44.csv cloud.ply --angles 0 --iteration 45

############################################################
## generate planar approximation
../globOptVis  --show3D  --pop-limit 3 -p primitives_it45.bonmin.csv -a points_primitives_it44.csv --title GlobOpt - [Dir-Colours] 45 iteration output --angle-gens 0 --draw-mode 1 --save-poly  --paral-colours --no-pts
