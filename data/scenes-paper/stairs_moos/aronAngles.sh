#!/bin/bash
rm ../bu_stairs_moos_angles/*.csv

scale=0.01;
recallThresh=$scale;
filterAmbig=1;
N=500000;
statLogPath="statLog.csv";

rm $statLogPath

/home/bontius/workspace/globOpt/globOpt/build/Release/bin/eval --stat-log $statLogPath --recall-thresh $recallThresh --n-rels $N --planes --assign gt.obj --cloud cloud.ply -p ../stairs_moos_globfit_0.5_300/primitives_it15.bonmin.csv -a ../stairs_moos_globfit_0.5_300/points_primitives_it14.csv --scale $scale --silent --filter-ambig $filterAmbig
cp ../stairs_moos_globfit_0.5_300/primitives_it15.bonmin.angles.csv ../bu_stairs_moos_angles/

/home/bontius/workspace/globOpt/globOpt/build/Release/bin/eval --stat-log $statLogPath --recall-thresh $recallThresh --n-rels $N --planes --assign gt.obj -p schnabel_minsup1500.primitives.csv -a schnabel_minsup1500.points_primitives.csv --cloud schnabel_minsup1500.cloud.ply --scale $scale --silent --filter-ambig $filterAmbig
cp schnabel_minsup1500.primitives.angles.csv ../bu_stairs_moos_angles/

#/home/bontius/workspace/globOpt/globOpt/build/Release/bin/eval --stat-log $statLogPath --recall-thresh $recallThresh --n-rels $N --planes --assign gt.obj --cloud cloud_100000.ply -p primitives.pearl.csv -a points_primitives.pearl.csv --scale $scale --silent --filter-ambig $filterAmbig
#cp primitives.pearl.angles.csv ../bu_stairs_moos_angles/

/home/bontius/workspace/globOpt/globOpt/build/Release/bin/eval --stat-log $statLogPath --recall-thresh $recallThresh --n-rels $N --planes --assign gt.obj --cloud cloud.ply --scale $scale -p globfit_0.5_300/patches.sub_0.5_300.csv -a globfit_0.5_300/points_patches.sub_0.5_300.csv --silent --filter-ambig $filterAmbig
cp patches.sub_0.5_300.angles.csv ../bu_stairs_moos_angles/

cp statLog.csv ../bu_stairs_moos_angles/
cat $statLogPath