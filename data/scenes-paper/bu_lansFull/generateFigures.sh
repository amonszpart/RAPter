#!/bin/sh

SCRIPT_PATH="../../../globOpt/scripts"

mkdir figure

############################################################
## generate stats from input


meshlabserver -i cloud.ply -o figure/cloud_pts.ply -s $SCRIPT_PATH/colorize.mlx -om vc vn

splatting figure/cloud_pts.ply figure/cloud.ply 1 0.0025


############################################################
## generate cleaned point-cloud model, and associated stats

meshlabserver -i figure/cloudRGBNormal.ply -o figure/cloud_cleaned_pts.ply -om vc vn

splatting figure/cloud_cleaned_pts.ply figure/cloud_cleaned.ply 1 0.0025

## planar approximation is not available yet
