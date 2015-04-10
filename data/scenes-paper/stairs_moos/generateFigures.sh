../compareToGlobfit.py -s 0.01 --primLimit 15 --run
../compareToGlobfit.py -s 0.01 --primLimit 300 --run

#!/bin/sh

SCRIPT_PATH="../../../../globOpt/scripts"

############################################################
## generate stats from input

python ../normal_distr.py cloud.ply ndistrCloudBinary.svg "stairs - Input" --noscatter

meshlabserver -i cloud.ply -o globfit_0.5_300/figure/cloud_pts.ply -s $SCRIPT_PATH/colorize.mlx -om vc vn

splatting globfit_0.5_300/figure/cloud_pts.ply globfit_0.5_300/figure/cloud.ply 1 0.005
cp globfit_0.5_300/figure/cloud_pts.ply globfit_0.5_15/figure/cloud_pts.ply
cp globfit_0.5_300/figure/cloud.ply globfit_0.5_15/figure/cloud.ply

############################################################
## globfit 15
../show.py -s 0.010000 -p globfit_0.5_15/patches.sub_0.5_15.csv -a globfit_0.5_15/points_patches.sub_0.5_15.csv --save-poly

cd globfit_0.5_15
mkdir figure

meshlabserver -i ../cloudRGBNormal_patches_reProj_noUnass.ply -o figure/globfit_15_pts_segments.ply -om vc vn
# input segments

splatting figure/globfit_15_pts_segments.ply figure/globfit_15_segments.ply 1 0.005

# output
meshlabserver -i cloudRGBNormal_patches_reProj_noUnass.ply -o figure/globfit_15_pts.ply -om vc vn
splatting figure/globfit_15_pts.ply figure/globfit_15.ply 1 0.005

#distribution
python ../normal_distr.py  cloudRGBNormal_patches_reProj_noUnass.ply ndistr_globfit_15.svg "Stairs - Globfit"

cp plane_mesh_patches.ply figure/globfit_15_planes.ply


############################################################
## globfit 300

../show.py -s 0.010000 -p globfit_0.5_300/patches.sub_0.5_300.csv -a globfit_0.5_300/points_patches.sub_0.5_300.csv --save-poly

cd globfit_0.5_300
mkdir figure

meshlabserver -i ../cloudRGBNormal_patches_reProj_noUnass.ply -o figure/globfit_300_pts_segments.ply -om vc vn
# input segments

splatting figure/globfit_300_pts_segments.ply figure/globfit_300_segments.ply 1 0.005

# output
meshlabserver -i cloudRGBNormal_patches_reProj_noUnass.ply -o figure/globfit_300_pts.ply -om vc vn
splatting figure/globfit_300_pts.ply figure/globfit_300.ply 1 0.005

#distribution
python ../normal_distr.py  cloudRGBNormal_patches_reProj_noUnass.ply ndistr_globfit_300.svg "Stairs - Globfit"

cp plane_mesh_patches.ply figure/globfit_300_planes.ply
