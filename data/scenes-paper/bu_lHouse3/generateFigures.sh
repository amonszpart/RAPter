#!/bin/sh

SCRIPT_PATH="../../../globOpt/scripts"

mkdir figure

############################################################
## generate stats from input

#distribution
python ../normal_distr.py cloud.ply ndistrCloudBinary.svg "House - Input" --noscatter

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





############################################################
## schnabel

../globOptVis --show3D --scale 0.02 --pop-limit 3 --title "Schnabel 10" --angle-gens 90 --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-scale --bg-colour .9,.9,.9 --no-rel -p schnabel_minsup10.primitives.csv -a schnabel_minsup10.points_primitives.csv --cloud schnabel_minsup10.cloud.ply  --perfect-angle 0.0001  --draw-mode 28 --save-poly 

## clean output name to avoid side effect with patches output
mv cloudRGBNormal_reProj_noUnass_noPrim.ply schnabel_minsup10_pts.ply

meshlabserver -i schnabel_minsup10_pts.ply -o figure/schnabel_minsup10_pts.ply -om vc vn

splatting figure/schnabel_minsup10_pts.ply figure/schnabel_minsup10.ply 2 0.015

#distribution
python ../normal_distr.py schnabel_minsup10_pts.ply ndistr_schnabel_minsup10.svg "lHouse - Schnabel (minsup=10)"

../globOptVis --show3D --scale 0.02 --pop-limit 3 --title "Schnabel 20000" --angle-gens 90 --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-scale --bg-colour .9,.9,.9 --no-rel -p schnabel_minsup10.primitives.csv -a schnabel_minsup10.points_primitives.csv --cloud schnabel_minsup10.cloud.ply  --perfect-angle 0.0001  --draw-mode 1 --save-poly  --no-pts 

mv plane_mesh.ply figure/schnabel_minsup10_planes.ply



############################################################
## pearl

../globOptVis --show3D --scale 0.02 --pop-limit 3 --title "Pearl" --angle-gens 90 --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-scale --bg-colour .9,.9,.9 --no-rel -p primitives.pearl.csv -a points_primitives.pearl.csv --perfect-angle 0.0001  --draw-mode 28 --save-poly 

## clean output name to avoid side effect with patches output
mv cloudRGBNormal_reProj_noUnass_noPrim.ply pearl_pts.ply

meshlabserver -i pearl_pts.ply -o figure/pearl_pts.ply -om vc vn

splatting figure/pearl_pts.ply figure/pearl.ply 2 0.015

#distribution
python ../normal_distr.py pearl_pts.ply ndistr_pearl.svg "lHouse - Pearl"

../globOptVis --show3D --scale 0.02 --pop-limit 3 --title "Pearl" --angle-gens 90 --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-scale --bg-colour .9,.9,.9 --no-rel -p primitives.pearl.csv -a points_primitives.pearl.csv --perfect-angle 0.0001  --draw-mode 1 --save-poly  --no-pts 

mv plane_mesh.ply figure/pearl_planes.ply
