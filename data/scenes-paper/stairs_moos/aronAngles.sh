#!/bin/bash
rm ../bu_stairs_moos_angles/*.csv

/home/bontius/workspace/globOpt/globOpt/build/Release/bin/eval --planes --assign gt.obj --cloud cloud.ply -p ../stairs_moos_globfit_0.5_300/primitives_it15.bonmin.csv -a ../stairs_moos_globfit_0.5_300/points_primitives_it14.csv --scale 0.01 --silent --filter-ambig 1
cp ../stairs_moos_globfit_0.5_300/primitives_it15.bonmin.angles.csv ../bu_stairs_moos_angles/

#/home/bontius/workspace/globOpt/globOpt/build/Release/bin/eval --planes --assign gt.obj -p schnabel_minsup10.primitives.csv -a schnabel_minsup10.points_primitives.csv --cloud schnabel_minsup10.cloud.ply --scale 0.01 --silent --filter-ambig 1
#cp schnabel_minsup10.primitives.angles.csv ../bu_stairs_moos_angles/

#/home/bontius/workspace/globOpt/globOpt/build/Release/bin/eval --planes --assign gt.obj --cloud cloud_100000.ply -p primitives.pearl.csv -a points_primitives.pearl.csv --scale 0.01 --silent --filter-ambig 1
#cp primitives.pearl.angles.csv ../bu_stairs_moos_angles/

/home/bontius/workspace/globOpt/globOpt/build/Release/bin/eval --planes --assign gt.obj --cloud cloud.ply --scale 0.010000 -p globfit_0.5_300/patches.sub_0.5_300.csv -a globfit_0.5_300/points_patches.sub_0.5_300.csv --silent --filter-ambig 1
cp globfit_0.5_300/patches.sub_0.5_300.angles.csv ../bu_stairs_moos_angles/