#!/bin/sh

SCRIPT_PATH="../../../globOpt/scripts"

mkdir figure

############################################################
## generate stats from input

python ../normal_distr.py cloud.ply ndistrCloudBinary.svg "Lans - Input"

meshlabserver -i cloud.ply -o figure/cloud_pts.ply -s $SCRIPT_PATH/colorize.mlx -om vc vn

splatting figure/cloud_pts.ply figure/cloud.ply 1 0.0025


############################################################
## generate stats for patches
../globOptVis  --show3D  --pop-limit 3 -p patches.csv -a points_primitives.csv --title GlobOpt - [Dir-Colours] patches --angle-gens 0 --draw-mode 28 --save-poly  --paral-colours 

## clean output name to avoid side effect with schnabel output
mv cloudRGBNormal_patches_reProj_noUnass_noPrim.ply patches.ply

#distribution
python ../normal_distr.py patches.ply ndistrPatches.svg "Lans - Patches"



############################################################
## generate cleaned point-cloud model, and associated stats

../globOptVis  --show3D  --pop-limit 3 -p primitives_it32.bonmin.csv -a points_primitives_it31.csv --title GlobOpt - [Dir-Colours] 32 iteration output --angle-gens 0 --draw-mode 28 --save-poly  --paral-colours 

#distribution
python ../normal_distr.py cloudRGBNormal_it32_reProj_noUnass_noPrim.ply ndistrIt32.svg "Lans - Iteration 32"

meshlabserver -i cloudRGBNormal_it32_reProj_noUnass_noPrim.ply -o figure/cloud_cleaned_pts.ply -om vc vn

splatting figure/cloud_cleaned_pts.ply figure/cloud_cleaned.ply 1 0.0025

#graph
python ../readGraphProperties.py primitives_it32.bonmin.csv points_primitives_it31.csv cloud.ply --angles 0 --iteration 32


############################################################
## generate planar approximation
../globOptVis  --show3D  --pop-limit 3 -p primitives_it32.bonmin.csv -a points_primitives_it31.csv --title GlobOpt - [Dir-Colours] 10 iteration output --angle-gens 0 --draw-mode 1 --save-poly  --paral-colours --no-pts 



############################################################
## schnabel

../globOptVis --show3D --scale 0.02 --pop-limit 3 --title "Schnabel 1000" --angle-gens 90 --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-scale --bg-colour .9,.9,.9 --no-rel -p schnabel_minsup1000.primitives.csv -a schnabel_minsup1000.points_primitives.csv --cloud schnabel_minsup1000.cloud.ply  --perfect-angle 0.0001  --draw-mode 28 --save-poly 

## clean output name to avoid side effect with patches output
mv cloudRGBNormal_reProj_noUnass_noPrim.ply schnabel_minsup1000_pts.ply

meshlabserver -i schnabel_minsup1000_pts.ply -o figure/schnabel_minsup1000_pts.ply -om vc vn

splatting figure/schnabel_minsup1000_pts.ply figure/schnabel_minsup1000.ply 1 0.0035

#distribution
python ../normal_distr.py schnabel_minsup1000_pts.ply ndistr_schnabel_minsup1000.svg "Lans - Schnabel (minsup=1000)"

../globOptVis --show3D --scale 0.02 --pop-limit 3 --title "Schnabel 1000" --angle-gens 90 --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-scale --bg-colour .9,.9,.9 --no-rel -p schnabel_minsup1000.primitives.csv -a schnabel_minsup1000.points_primitives.csv --cloud schnabel_minsup1000.cloud.ply  --perfect-angle 0.0001  --draw-mode 1 --save-poly  --no-pts 

mv plane_mesh.ply figures/schnabel_minsup1000_planes.ply



############################################################
## pearl 1

../globOptVis --show3D --scale 0.02 --pop-limit 3 --title "Pearl1" --angle-gens 90 --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-scale --bg-colour .9,.9,.9 --no-rel -p primitives.pearl.1.csv -a points_primitives.pearl.1.csv --perfect-angle 0.0001  --draw-mode 28 --save-poly  --cloud cloud.pearl.100000.ply

## clean output name to avoid side effect with patches output
mv cloudRGBNormal_reProj_noUnass_noPrim.ply pearl1_pts.ply

meshlabserver -i pearl1_pts.ply -o figure/pearl1_pts.ply -om vc vn

splatting figure/pearl1_pts.ply figure/pearl1.ply 1 0.002

#distribution
python ../normal_distr.py pearl1_pts.ply ndistr_pearl1.svg "Nola - Pearl"

../globOptVis --show3D --scale 0.02 --pop-limit 3 --title "Pearl1" --angle-gens 90 --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-scale --bg-colour .9,.9,.9 --no-rel -p primitives.pearl.1.csv -a points_primitives.pearl.1.csv --perfect-angle 0.0001  --draw-mode 1 --save-poly  --no-pts --cloud cloud.pearl.100000.ply

mv plane_mesh.ply figure/pearl1_planes.ply

############################################################
## globfit
../compareToGlobfit.py -s 0.001 --primLimit 50 --run
../compareToGlobfit.py -s 0.001 --primLimit 100 --run
../compareToGlobfit.py -s 0.001 --primLimit 150 --run

