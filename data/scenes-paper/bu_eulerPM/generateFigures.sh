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

#meshlabserver -i cloudRGBNormal_patches_reProj_noUnass_noPrim.ply -o figure/cloud_cleaned.ply -om vc vn

splatting figure/cloud_cleaned.ply figure/cloud_cleaned.ply 0 0.001


############################################################
## generate cleaned point-cloud model, and associated stats
../globOptVis  --show3D  --pop-limit 3 -p primitives_it69.bonmin.csv -a points_primitives_it68.csv --title GlobOpt - [Dir-Colours] 69 iteration output --angle-gens 0 --draw-mode 28 --save-poly  --paral-colours 

#distribution
python ../normal_distr.py cloudRGBNormal_it69_reProj_noUnass_noPrim.ply ndistrIt69.svg "Euler - Iteration 69"

#graph
python ../readGraphProperties.py primitives_it69.bonmin.csv points_primitives_it68.csv cloud.ply --angles 0 --iteration 69

meshlabserver -i cloudRGBNormal_it69_reProj_noUnass_noPrim.ply -o figure/cloud_cleaned_pts.ply -om vc vn

splatting figure/cloud_cleaned_pts.ply figure/cloud_cleaned.ply 0 0.001

############################################################
## generate planar approximation
../globOptVis  --show3D  --pop-limit 3 -p primitives_it69.bonmin.csv -a points_primitives_it68.csv --title GlobOpt - [Dir-Colours] 69 iteration output --angle-gens 0 --draw-mode 1 --save-poly  --paral-colours --no-pts --no-rel --no-clusters --no-pop --no-scale  --no-ids








############################################################
## schnabel

../globOptVis --show3D --scale 0.02 --pop-limit 3 --title "Schnabel 50" --angle-gens 90 --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-scale --bg-colour .9,.9,.9 --no-rel -p schnabel_minsup50.primitives.csv -a schnabel_minsup50.points_primitives.csv --cloud schnabel_minsup50.cloud.binary.ply  --perfect-angle 0.0001  --draw-mode 28 --save-poly 

## clean output name to avoid side effect with patches output
mv cloudRGBNormal_reProj_noUnass_noPrim.ply schnabel_minsup50_pts.ply

meshlabserver -i schnabel_minsup50_pts.ply -o figure/schnabel_minsup50_pts.ply -om vc vn

splatting figure/schnabel_minsup50_pts.ply figure/schnabel_minsup50.ply 0 0.001

#distribution
python ../normal_distr.py schnabel_minsup50_pts.ply ndistr_schnabel_minsup50.svg "eulerPM - Schnabel (minsup=50)"

python ../normal_distr.py schnabel_minsup50_pts.ply ndistr_schnabel_minsup50_log.svg "eulerPM - Schnabel (minsup=50)" --applylog

../globOptVis --show3D --scale 0.02 --pop-limit 3 --title "Schnabel 50" --angle-gens 90 --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-scale --bg-colour .9,.9,.9 --no-rel -p schnabel_minsup50.primitives.csv -a schnabel_minsup50.points_primitives.csv --cloud schnabel_minsup50.cloud.binary.ply  --perfect-angle 0.0001  --draw-mode 1 --save-poly  --no-pts 

mv plane_mesh.ply figure/schnabel_minsup50_planes.ply

