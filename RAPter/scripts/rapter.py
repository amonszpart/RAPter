#!/usr/bin/python

import argparse
import os
import sys          # exit
import shutil       # copyfile
import math         # pi
import subprocess   # Popen


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
        ret = os.system(cmd) >> 8 # python thing
        if ret != 0:
            if not noExit:
                print("call returned error ", ret, ", aborting")
                sys.exit()
        return ret

parser = argparse.ArgumentParser()

parser.add_argument( "-s"  , "--scale"      , dest="scale"        , type=float, default=0.05, help="Scale (rho) parameter, the smallest feature size to preserve [0.001..0.05]")
parser.add_argument( "--al", "--angle-limit", dest="angleLimit"   , type=float, default=15  , help="Angle threshlod (tau) parameter in degrees [5..45]")
parser.add_argument( "--pw", "--pairwise"   , dest="pw"           , type=float, default=1.0 , help="Weight of pairwise term [0.1..10^6]" )
parser.add_argument( "--ag", "--angle-gens" , dest="angleGens"    , type=float, default=[0,90], help="Weight of pairwise term [0.1..10^6]", action="append" )
parser.add_argument( "--it", "--iterations" , dest="nbExtraIterations", type=int, default=15, help="How many iterations to run [5..20]")

parser.add_argument( "--pl", "--popLimit"   , dest="popLimit"     , type=int  , default=5   , help="Filters primitives having less than this many points assigned [3..100]")
parser.add_argument( "--sp", "--spatial"    , dest="spatial"      , type=float,               help="Weight of spatial term [0.1, pw/10., pw/5., pw/2.]" )
parser.add_argument( "-l"  , "--lines"      , dest="lines"        , action="store_true"     , help="Work in 2D with lines instead of planes." )
parser.add_argument("--vl" , "--var-limit"  , dest="variableLimit", type=int  , default=1000, help="Maximum number of variables (primitives) for the optimisation. [500..3000]")

parser.add_argument( "-d"  , "--data"       , dest="data"         , type=float, default=1e5 , help="Weight of data term [10^5, 10^6]" )
parser.add_argument( "-p"  , "--primitives" , dest="primitives"   , type=str  ,               help="Input primitives, e.g. existing segmentation segments.csv" )
parser.add_argument( "-a"  , "--assoc"      , dest="associations" , type=str  ,               help="Input point-primitive associations, e.g. existing segmentation's points_segments.csv" )

parser.add_argument( "--dry", action="store_true"                 , help="Call the scripts (disabled by default)" )
parser.add_argument( "--vis", action="store_false", default = True, help="Enable visualization" )

parser.add_argument("--segment-scale-mult", dest="segmentScaleMultiplier", type=float, default=1.0, help="Multiply scale by this value for the segmentation step. [0.5, 1.0, 2.0]")
parser.add_argument("--ald", "--angle-limit-divisor", dest="angleLimitDivisor", type=float, default=1.0, help="Divide angle threshold (tau) by this number for candidate generation. [2.0, 1.0, 0.5]")

args = parser.parse_args()

if not os.path.isfile("cloud.ply"):
    print("Need \"cloud.ply\" to exist, assuming it's the pointcloud");
    sys.exit(1);

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
    setattr( args, "flag3D"     ,"3D"             );
    setattr( args, "tripletSafe","--triplet-safe" ); # Values: ["", "--triplet-safe"]
else:
    setattr( args, "flag3D"     ,"" );
    setattr( args, "tripletSafe","" );

########################################################################################################################
# Do segmentation
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

########################################################################################################################

########################################################################################################################

angleGens           = "0"
primitives          = "patches.csv"
associations        = "points_primitives.csv"
keepSingles         = "--keep-singles"
allowPromoted       = "--allow-promoted"
smallThreshMult     = 2.
promRem             = 0                         # remaining primitives to promote
collapseThreshDeg   = 0.4                       # initialize optimisation with the closest two orientations merged, if their difference is < colleapseThreshDeg degrees.
algCode             = 0                         # 0==B_BB, OA, QG, 3==Hyb, ECP, IFP
adopt               = 0

iteration = 0

# (2) Generate - generate candidates
cmd = "%s --generate%s -sc %f -al %f -ald %f --patch-pop-limit %d -p %s --assoc %s --angle-gens %s --small-thresh-mult %f --var-limit %d --small-mode 0 %s %s %s" \
        % (rapterExec, args.flag3D, args.scale, args.angleLimit, args.angleLimitDivisor, args.popLimit, primitives, associations, angleGens, \
           smallThreshMult, args.variableLimit, \
           args.tripletSafe, keepSingles, allowPromoted );
promRem = call( cmd, args.dry, True );
print "[rapter.py] Remaining smalls to promote: ", promRem
# Don't decrase area threshold until no more candidates to promote
if promRem != 0:
    decreaseLevel = False; # set to true each iteration

# (3) Formulate - create optimisation problem
cmd = "%s --formulate%s --scale %f --unary %f --pw %f --spat-weight %f --spat-dist-mult 2. --patch-pop-limit %d --angle-gens %s --cloud cloud.ply --candidates candidates_it%d.csv -a %s --collapse-angle-deg %f --trunc-angle %f --constr-mode patch --dir-bias 0 --no-clusters --cmp 0 --freq-weight 0 --cost-fn spatsqrt" \
       % ( rapterExec, args.flag3D, args.scale, args.data, args.pw, args.spatial, args.popLimit, angleGens, iteration, associations, collapseThreshDeg, args.angleLimit)
call( cmd, args.dry );

# (4) Solve
cmd = "%s --solver%s bonmin --problem problem -v --time -1 --bmode %d --angle-gens %s --candidates candidates_it%d.csv" \
       % ( rapterExec, args.flag3D, algCode, angleGens, iteration );
call( cmd, args.dry )

if args.vis:
    show( "primitives_it%d.bonmin.csv" % iteration, associations, "\"RAPter - Iteration%d\"" % iteration, args );

# (6) CoPlanarity
cmd = "%s --merge%s --scale %f --adopt %d --prims primitives_it%d.bonmin.csv -a %s --angle-gens %s --patch-pop-limit %d" \
      % ( rapterExec, args.flag3D, args.scale, adopt, iteration, associations, angleGens, args.popLimit );
call( cmd, args.dry )

# Don't copy promoted patches' directions to other patches after 4 iterations (c==3), since they are not reliable anymore
if iteration == 3:
    allowPromoted = ""

# Don't throw away single directions before the 3rd (c==1) iteration. 
# This will keep large patches, even if they don't copy to anywhere for later.
if iteration == 1:
    keepSingles="";

# If we are still promoting small patches on this working scale, make sure to run more iterations
if iteration == args.nbExtraIterations and promRem != 0:
    args.nbExtraIterations += 1;

# Increment iteration counter (while loop)
iteration += 1;
