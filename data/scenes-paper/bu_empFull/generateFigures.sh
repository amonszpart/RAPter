#!/bin/sh

SCRIPT_PATH="../../../globOpt/scripts"

ln -s cloudBinary.ply cloud.ply
mkdir figure

############################################################
## generate stats from input

#distribution
python ../normal_distr.py cloud.ply ndistrCloudBinary.svg "Empire - Input"

meshlabserver -i cloud.ply -o figure/cloud_pts.ply -s $SCRIPT_PATH/colorize.mlx -om vc vn

splatting figure/cloud_pts.ply figure/cloud.ply 1 0.001


############################################################
## generate stats for patches
../globOptVis  --show3D  --pop-limit 3 -p patches.csv -a points_primitives.csv --title GlobOpt - [Dir-Colours] patches --angle-gens 0 --draw-mode 28 --save-poly  --paral-colours 

#distribution
python ../normal_distr.py cloudRGBNormal_patches_reProj_noUnass_noPrim.ply ndistrPatches.svg "Empire - Patches"


############################################################
## generate cleaned point-cloud model, and associated stats
../globOptVis  --show3D  --pop-limit 3 -p primitives_it10.bonmin.csv -a points_primitives_it9.csv --title GlobOpt - [Dir-Colours] 10 iteration output --angle-gens 0 --draw-mode 28 --save-poly  --paral-colours 

#distribution
python ../normal_distr.py cloudRGBNormal_it10_reProj_noUnass_noPrim.ply ndistrIt10.svg "Empire - Iteration 10"

meshlabserver -i cloudRGBNormal_it10_reProj_noUnass_noPrim.ply -o figure/cloud_cleaned_pts.ply -om vc vn

splatting figure/cloud_cleaned_pts.ply figure/cloud_cleaned.ply 1 0.002

#graph
python ../readGraphProperties.py primitives_it10.bonmin.csv points_primitives_it9.csv cloud.ply --angles 0 --iteration 10

############################################################
## generate planar approximation
../globOptVis  --show3D  --pop-limit 3 -p primitives_it10.bonmin.csv -a points_primitives_it9.csv --title GlobOpt - [Dir-Colours] 10 iteration output --angle-gens 0 --draw-mode 1 --save-poly  --paral-colours --no-pts






############################################################
## schnabel

../globOptVis --show3D --scale 0.02 --pop-limit 3 --title "Schnabel 1000" --angle-gens 90 --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-scale --bg-colour .9,.9,.9 --no-rel -p schnabel_minsup20000.primitives.csv -a schnabel_minsup20000.points_primitives.csv --cloud schnabel_minsup20000.cloud.ply  --perfect-angle 0.0001  --draw-mode 28 --save-poly 

## clean output name to avoid side effect with patches output
mv cloudRGBNormal_reProj_noUnass_noPrim.ply schnabel_minsup20000_pts.ply

meshlabserver -i schnabel_minsup20000_pts.ply -o figure/schnabel_minsup20000_pts.ply -om vc vn

splatting figure/schnabel_minsup20000_pts.ply figure/schnabel_minsup20000.ply 1 0.0035

#distribution
python ../normal_distr.py schnabel_minsup20000_pts.ply ndistr_schnabel_minsup20000.svg "Lans - Schnabel (minsup=20000)"

../globOptVis --show3D --scale 0.02 --pop-limit 3 --title "Schnabel 20000" --angle-gens 90 --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-scale --bg-colour .9,.9,.9 --no-rel -p schnabel_minsup20000.primitives.csv -a schnabel_minsup20000.points_primitives.csv --cloud schnabel_minsup20000.cloud.ply  --perfect-angle 0.0001  --draw-mode 1 --save-poly  --no-pts 

mv plane_mesh.ply figure/schnabel_minsup20000_planes.ply




############################################################
## pearl

../globOptVis --show3D --scale 0.02 --pop-limit 3 --title "Pearl" --angle-gens 90 --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-scale --bg-colour .9,.9,.9 --no-rel -p primitives.pearl.csv -a points_primitives.pearl.csv --perfect-angle 0.0001  --draw-mode 28 --save-poly  --cloud cloud_100000.ply 

## clean output name to avoid side effect with patches output
mv cloudRGBNormal_reProj_noUnass_noPrim.ply pearl_pts.ply

meshlabserver -i pearl_pts.ply -o figure/pearl_pts.ply -om vc vn

splatting figure/pearl_pts.ply figure/pearl.ply 1 0.0035

#distribution
python ../normal_distr.py pearl_pts.ply ndistr_pearl.svg "lHouse - Pearl"

../globOptVis --show3D --scale 0.02 --pop-limit 3 --title "Pearl" --angle-gens 90 --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-scale --bg-colour .9,.9,.9 --no-rel -p primitives.pearl.csv -a points_primitives.pearl.csv --perfect-angle 0.0001  --draw-mode 1 --save-poly  --no-pts   --cloud cloud_100000.ply 

mv plane_mesh.ply figure/pearl_planes.ply


############################################################
## globfit 13

../show.py -s 0.02 -p globfit_0.5_13/patches.sub_0.5_13.csv -a globfit_0.5_13/points_patches.sub_0.5_13.csv --save-poly

cd globfit_0.5_13
mkdir figure

meshlabserver -i ../cloudRGBNormal_patches_reProj_noUnass.ply -o figure/globfit_13_pts_segments.ply -om vc vn
# input segments

splatting figure/globfit_13_pts_segments.ply figure/globfit_13_segments.ply 1 0.003

# output
meshlabserver -i cloudRGBNormal_patches_reProj_noUnass.ply -o figure/globfit_13_pts.ply -om vc vn
splatting figure/globfit_13_pts.ply figure/globfit_13.ply 1 0.003


############################################################
## lafarge

../show.py -s 0.02 -p primitives.lafarge.csv -a points_primitives.lafarge.csv  --save-poly

# output
meshlabserver -i cloudRGBNormal_reProj_noUnass.ply -o figure/lafarge_pts.ply -om vc vn
splatting figure/lafarge_pts.ply figure/lafarge.ply 1 0.003

#distribution
python ../normal_distr.py figure/lafarge_pts.ply ndistr_lafarge.svg "empire - Lafarge"

cp plane_mesh_lafarge.ply figure/lafarge_planes.ply
