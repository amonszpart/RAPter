#!/bin/bash
# Run script from scenes or scenes-paper to create all symlinks

rapterRoot="/home/bontius/workspace/globOpt"

ln -sf $rapterRoot/RAPter/build/Release/bin/rapter
ln -sf $rapterRoot/RAPter/build/Debug/bin/rapter rapterd

ln -sf $rapterRoot/visualization/build/Release/bin/rapterVis
ln -sf $rapterRoot/visualization/build/Debug/bin/rapterVis raperVisd

ln -sf $rapterRoot/RAPter/build/Release/bin/corresp
ln -sf $rapterRoot/RAPter/scripts/divide.py

ln -sf $rapterRoot/RAPter/scripts/run.sh
ln -sf $rapterRoot/RAPter/scripts/runRepr.sh

ln -sf $rapterRoot/RAPter/build/Release/bin/pearl
ln -sf $rapterRoot/RAPter/build/Debug/bin/pearl pearld
ln -sf $rapterRoot/RAPter/build/Release/bin/ransac
ln -sf $rapterRoot/RAPter/build/Debug/bin/ransac ransacd
ln -sf $rapterRoot/RAPter/build/Release/bin/toGlobFit

ln -sf $rapterRoot/RAPter/build/Release/bin/refit
ln -sf $rapterRoot/RAPter/scripts/refit.py

ln -sf $rapterRoot/RAPter/scripts/show.py
ln -sf $rapterRoot/RAPter/scripts/runGlobfit.py
ln -sf $rapterRoot/RAPter/scripts/compareToGlobfit.py
ln -sf $rapterRoot/RAPter/scripts/noisePw2Html.py
ln -sf $rapterRoot/RAPter/scripts/noisePwMatrix.py
ln -sf $rapterRoot/RAPter/scripts/runSegmentation.py

ln -sf $rapterRoot/evaluation/normal_distr.py
ln -sf $rapterRoot/evaluation/readGraphProperties.py
ln -sf $rapterRoot/evaluation/collectStatistics.py
ln -sf $rapterRoot/evaluation/computeNormalEstimationError.py
ln -sf $rapterRoot/evaluation/generatePolarPrimDirections.py

ln -sf $rapterRoot/RAPter/scripts/showcands.sh
ln -sf $rapterRoot/RAPter/scripts/showPearl.sh
ln -sf $rapterRoot/RAPter/scripts/showquick.sh
ln -sf $rapterRoot/RAPter/scripts/cleanupFolder.sh