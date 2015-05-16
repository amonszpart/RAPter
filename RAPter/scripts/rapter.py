#!/usr/bin/python

import argparse
import os
import shutil # copyfile
import math # pi
import subprocess


rapterRoot = "/home/bontius/workspace/globOpt";
rapterExec = os.path.join( rapterRoot, "RAPter", "build", "Release", "bin", "rapter" );

def show( primitivesPath, associationsPath, title, args ):
    cmd = os.path.join("..","rapterVis --show%s --scale %f --pop-limit %d -p %s -a %s --title %s --angle-gens %s --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-rel --no-scale --bg-colour 1.,1.,1. --no-rel" \
            % ( args.flag3D, args.scale, args.popLimit, primitivesPath, associationsPath, title, args.angleGensStr ) );
    print cmd
    #print os.spawnlp( os.P_NOWAIT, "..", cmd )
    subprocess.Popen( cmd, shell=True );


def call( cmd, dry = True, noExit = False ):
    print("%s" % (cmd))
    if dry:
        print "DRY"
    else:
        print "RUN"
    if not dry:
        ret = os.system(cmd)
        if ret != 0:
            print("call returned error ", ret, ", aborting")
            if not noExit:
                exit
            else:
                print("call returned ", ret)

parser = argparse.ArgumentParser()

parser.add_argument( "-s"  , "--scale"      , dest="scale"       , type=float, default=0.05, help="Scale (rho) parameter, the smallest feature size to preserve [0.001..0.05]")
parser.add_argument( "--al", "--angle-limit", dest="angleLimit"  , type=float, default=15  , help="Angle threshlod (tau) parameter in degrees [5..45]")
parser.add_argument( "--pw", "--pairwise"   , dest="pw"          , type=float, default=1.0 , help="Weight of pairwise term [0.1..10^6]" )
parser.add_argument( "--ag", "--angle-gens" , dest="angleGens"   , type=float, default=[0,90], help="Weight of pairwise term [0.1..10^6]", action="append" )

parser.add_argument( "--pl", "--popLimit"   , dest="popLimit"    , type=int  , default=5   , help="Filters primitives having less than this many points assigned [3..100]")
parser.add_argument( "--sp", "--spatial"    , dest="spatial"     , type=float,               help="Weight of spatial term [0.1, pw/10., pw/5., pw/2.]" )
parser.add_argument( "-l"  , "--lines"      , dest="lines"       , action="store_true"     , help="Work in 2D with lines instead of planes." )

parser.add_argument( "-d"  , "--data"       , dest="data"        , type=float, default=1e5 , help="Weight of data term [10^5, 10^6]" )
parser.add_argument( "-p"  , "--primitives" , dest="primitives"  , type=str  ,               help="Input primitives, e.g. existing segmentation segments.csv" )
parser.add_argument( "-a"  , "--assoc"      , dest="associations", type=str  ,               help="Input point-primitive associations, e.g. existing segmentation's points_segments.csv" )

parser.add_argument( "--dry", action="store_true"                 , help="Call the scripts (disabled by default)" )
parser.add_argument( "--vis", action="store_false", default = True, help="Enable visualization" )

parser.add_argument("--segment-scale-mult", dest="segmentScaleMultiplier", type=float, default=1.0, help="Multiply scale by this value for the segmentation step. [0.5, 1.0, 2.0]")

args = parser.parse_args()

# if not args.scale:
#     print("Need scale -s, --scale")
#     exit
# if not args.angleLimit:
#     print("Need angleLimit! Set using '-al' or '--angle-limit'!")
#     exit

# convert to radians
args.angleLimit = args.angleLimit / 180.0 * math.pi
args.angleGensStr = ",".join( str(e) for e in args.angleGens )

print( "--popLimit %d \twill keep all primitives, that have more than this number of assigned points" % (args.popLimit) );
if not args.spatial:
    args.spatial = args.pw / 10.
    print( "--spatial %.3f" % (args.spatial) );

if not args.lines:
    setattr( args, "flag3D","3D");
else:
    setattr(args,"flag3D","");

# (0) Do segmentation
if not args.primitives or not args.associaitons:
    cmd = "%s --segment%s --scale %f --angle-limit %f --angle-gens %s --patch-pop-limit %d --dist-limit-mult %f" \
          % ( rapterExec, args.flag3D, args.scale, args.angleLimit, args.angleGensStr, args.popLimit, args.segmentScaleMultiplier )
    call( cmd, args.dry );
    # # save output
    if ( os.path.isfile("patches.csv") and os.path.isfile("points_primitives.csv") ):
        if os.path.isfile("segments.csv"):
            shutil.copyfile( "segments.csv", "segments.csv.bak" );
        shutil.copyfile( "patches.csv", "segments.csv" )
        if os.path.isfile("points_segments.csv"):
            shutil.copyfile( "points_segments.csv", "points_segments.csv.bak" );
        shutil.copyfile( "points_primitives.csv", "points_segments.csv" )

    if args.vis:
        show( "segments.csv", "points_segments.csv", "\"RAPter - Segmentation\"", args );