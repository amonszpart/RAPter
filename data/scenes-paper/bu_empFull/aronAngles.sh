#!/bin/bash
/home/bontius/workspace/globOpt/globOpt/build/Release/bin/eval --planes --assign gt.obj --cloud cloud.ply -p primitives_it11.bonmin.csv -a points_primitives_it10.csv --scale 0.0025 --silent

/home/bontius/workspace/globOpt/globOpt/build/Release/bin/eval --planes --assign gt.obj --cloud schnabel_minsup20000.cloud.ply -p schnabel_minsup20000.primitives.csv -a schnabel_minsup20000.points_primitives.csv --scale 0.0025 --silent

/home/bontius/workspace/globOpt/globOpt/build/Release/bin/eval --planes --assign gt.obj --cloud cloud_100000.ply -p primitives.pearl.csv -a points_primitives.pearl.csv --scale 0.0025 --silent

/home/bontius/workspace/globOpt/globOpt/build/Release/bin/eval --planes --assign gt.obj --cloud cloud.ply -p globfit_0.5_13/patches.sub_0.5_13.csv -a globfit_0.5_13/points_patches.sub_0.5_13.csv --scale 0.0025 --silent