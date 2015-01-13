#!/bin/sh

ln -s cloudBinary.ply cloud.ply

############################################################
## generate stats from input

#distribution
python ../normal_distr.py cloudBinary.ply ndistrCloudBinary.svg "Euler - Input"

############################################################
## generate stats for patches
../globOptVis  --show3D  --pop-limit 3 -p patches.csv -a points_primitives.csv --title GlobOpt - [Dir-Colours] patches --angle-gens 0 --draw-mode 28 --save-poly  --paral-colours 

#distribution
python ../normal_distr.py cloudRGBNormal_patches_reProj_noUnass_noPrim.ply ndistrPatches.svg "Euler - Patches"

############################################################
## generate cleaned point-cloud model, and associated stats
../globOptVis  --show3D  --pop-limit 3 -p primitives_it69.bonmin.csv -a points_primitives_it68.csv --title GlobOpt - [Dir-Colours] 69 iteration output --angle-gens 0 --draw-mode 28 --save-poly  --paral-colours 

#distribution
python ../normal_distr.py cloudRGBNormal_it69_reProj_noUnass_noPrim.ply ndistrIt69.svg "Euler - Iteration 69"

#graph
python ../readGraphProperties.py primitives_it69.bonmin.csv points_primitives_it68.csv cloud.ply --angles 0 --iteration 69

############################################################
## generate planar approximation
../globOptVis  --show3D  --pop-limit 3 -p primitives_it69.bonmin.csv -a points_primitives_it68.csv --title GlobOpt - [Dir-Colours] 69 iteration output --angle-gens 0 --draw-mode 28 --save-poly  --paral-colours --no-pts
