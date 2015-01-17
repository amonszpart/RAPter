#!/bin/sh

SCRIPT_PATH="../../../globOpt/scripts"

ln -s cloudBinary.ply cloud.plymkdir figure

############################################################
## generate stats from input

#distribution
python ../normal_distr.py cloud.ply ndistrCloudBinary.svg "Nola - Input"

meshlabserver -i cloud.ply -o figure/cloud_pts.ply -s $SCRIPT_PATH/colorize.mlx -om vc vn

splatting figure/cloud_pts.ply figure/cloud.ply 1 0.002


############################################################
## generate stats for patches
../globOptVis  --show3D  --pop-limit 3 -p patches.csv -a points_primitives.csv --title GlobOpt - [Dir-Colours] patches --angle-gens 0 --draw-mode 28 --save-poly  --paral-colours 

#distribution
python ../normal_distr.py cloudRGBNormal_patches_reProj_noUnass_noPrim.ply ndistrPatches.svg "Nola - Patches"


############################################################
## generate cleaned point-cloud model, and associated stats
../globOptVis  --show3D  --pop-limit 3 -p primitives_it26.bonmin.csv -a points_primitives_it25.csv --title GlobOpt - [Dir-Colours] 10 iteration output --angle-gens 0 --draw-mode 28 --save-poly  --paral-colours 

#distribution
python ../normal_distr.py cloudRGBNormal_it26_reProj_noUnass_noPrim.ply ndistrIt26.svg "Nola - Iteration 26"

meshlabserver -i cloudRGBNormal_it26_reProj_noUnass_noPrim.ply -o figure/cloud_cleaned_pts.ply -om vc vn

splatting figure/cloud_cleaned_pts.ply figure/cloud_cleaned.ply 1 0.002

#graph
python ../readGraphProperties.py primitives_it26.bonmin.csv points_primitives_it25.csv cloud.ply --angles 0 --iteration 10

############################################################
## generate planar approximation
../globOptVis  --show3D  --pop-limit 3 -p primitives_it26.bonmin.csv -a points_primitives_it25.csv --title GlobOpt - [Dir-Colours] 25 iteration output --angle-gens 0 --draw-mode 1 --save-poly  --paral-colours --no-pts
