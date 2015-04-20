#!/bin/bash
rm ../bu_lHouse3D_angles/*.csv

/home/bontius/workspace/globOpt/globOpt/build/Release/bin/eval --planes --assign gt.obj --cloud cloud.ply -p primitives_it20.bonmin.csv -a points_primitives_it19.csv --scale 0.0175 --silent --filter-ambig 1
cp primitives_it20.bonmin.angles.csv ../bu_lHouse3D_angles/

/home/bontius/workspace/globOpt/globOpt/build/Release/bin/eval --planes --assign gt.obj -p schnabel_minsup10.primitives.csv -a schnabel_minsup10.points_primitives.csv --cloud schnabel_minsup10.cloud.ply --scale 0.0175 --silent --filter-ambig 1
cp schnabel_minsup10.primitives.angles.csv ../bu_lHouse3D_angles/
wc -l schnabel_minsup10.primitives.angles.csv

/home/bontius/workspace/globOpt/globOpt/build/Release/bin/eval --planes --assign gt.obj --cloud cloud.ply -p primitives.pearl.csv -a points_primitives.pearl.csv --scale 0.0175 --silent --filter-ambig 1
cp primitives.pearl.angles.csv ../bu_lHouse3D_angles/
wc -l primitives.pearl.angles.csv

/home/bontius/workspace/globOpt/globOpt/build/Release/bin/eval --planes --assign gt.obj --cloud ../bu_lHouse3D_globfit/cloud.ply -p ../bu_lHouse3D_globfit/primitives.globfit.csv -a ../bu_lHouse3D_globfit/points_primitives.globfit.csv --scale 0.0175 --silent --filter-ambig 1
cp ../bu_lHouse3D_globfit/primitives.globfit.angles.csv ../bu_lHouse3D_angles/