#!/usr/bin/python

import argparse
import os
import sys          # exit
import shutil       # copyfile
import math         # pi
import subprocess   # Popen


rapterRoot = "/home/bontius/workspace/RAPter/";
rapterExec = os.path.join( rapterRoot, "RAPter", "build", "Release", "bin", "rapter" );

def show( primitivesPath, associationsPath, title, args ):
    cmd = os.path.join("..","rapterVis --show%s --scale %f --pop-limit %d -p %s -a %s --cloud %s --title %s --angle-gens %s --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-rel --no-scale --bg-colour 1.,1.,1. --no-rel" \
            % ( args.flag3D, args.scale, args.popLimit, primitivesPath, associationsPath, args.cloud, title, args.angleGensStr ) );
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

def runRepr( rprPrims, rprAssoc, rprIter, args, angleGens, keepSingles ):
    print( "Running with prims: %s, assoc: %s, iteration: %d" % (rprPrims,rprAssoc,rprIter) );

    rprRepr      = "representatives_it%d.csv" % rprIter             # representatives output, contains one patch for each dId
    rprReprAssoc = "points_representatives_it%d.csv" % rprIter      # representatives output, contains associations for representative primitives only
    rprCands     = "candidates_representatives_it%d.csv" % rprIter  # candidates generated from representatives
    rprReprOpt   = "representatives_it%d.bonmin.csv" % rprIter      # new representatives chosen from candidates
    rprPrimBak   =  "%s.lvl1.csv" % os.path.splitext(rprPrims)[0]   # "`cutExt $rprPrims`".lvl1.csv


    rprNextId    = rprIter + 1; #`expr $c + 1`; #candidates will output here automatically...so we need to know
    rprPw        = args.pw
    rprAngLimit  = args.angleLimit

    # representatives
    cmd = "%s --represent%s -p %s -a %s -sc %f --cloud %s --angle-gens %s" \
           % (args.rapterExec, args.flag3D, rprPrims, rprAssoc, args.scale, args.cloud, angleGens );
    #my_exec "$executable --represent$flag3D -p $rprPrims -a $rprAssoc -sc $scale --cloud cloud.ply --angle-gens $anglegens"
    call( cmd, args.dry );

    if not args.dry:
        #echo "mv representatives.csv $rprRepr"
        #mv representatives.csv $rprRepr
        shutil.move("representatives.csv",rprRepr);
        #echo "mv points_representatives.csv $rprReprAssoc"
        #mv points_representatives.csv $rprReprAssoc
        shutil.move("points_representatives.csv",rprReprAssoc);


    # ShowRepr
    #cmd="../globOptVis --show$flag3D--scale $scale --pop-limit $poplimit --title \"Representatives\" --angle-gens $angleGens --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-scale --bg-colour .9,.9,.9 --ids --no-rel -p $repr -a $assoc $"
    #my_exec "../globOptVis --show$flag3D --scale $scale --pop-limit $poplimit -p $rprRepr -a $rprReprAssoc --title \"Representatives\" $visdefparam &"

    # Generate from Repr
    if not args.dry:
        #echo "mv candidates_it${rprNextId}.csv candidates_it${rprNextId}_tmp.csv" # move tmp out of the way
        #mv candidates_it${rprNextId}.csv candidates_it${rprNextId}_tmp.csv # move tmp out of the way
        if os.path.isfile( "candidates_it%d.csv" % rprNextId ):
            shutil.move( "candidates_it%d.csv" % rprNextId, "candidates_it%d_tmp.csv" % rprNextId );

    #cmd = "$executable --generate$flag3D $tripletSafe -sc $scale -al $rprAngLimit -ald ${cand_anglediv} --small-mode 0 --patch-pop-limit $poplimit --angle-gens $candAngleGens --small-thresh-mult $smallThresh -p $rprRepr --assoc $rprReprAssoc --keep-singles"
    #my_exec "%s --generate%s -sc $scale -al $rprAngLimit -ald 1.0 --patch-pop-limit $poplimit -p $rprRepr --assoc $rprReprAssoc --angle-gens $candAngleGens --small-thresh-mult %f --small-mode 0 %s %s"
    cmd = "%s --generate%s -sc %f --cloud %s -al %f -ald 1.0 --patch-pop-limit %d -p %s --assoc %s --angle-gens %s --small-thresh-mult %f --small-mode 0 %s %s" \
        % ( args.rapterExec, args.flag3D, args.scale, args.cloud, rprAngLimit, args.popLimit, rprRepr, rprReprAssoc, candAngleGens, \
           args.smallThreshMult, args.tripletSafe, keepSingles );
    call( cmd, args.dry );

    if not args.dry:
        # echo "mv candidates_it${rprNextId}.csv $rprCands"
        # mv candidates_it${rprNextId}.csv $rprCands
        shutil.move( "candidates_it%d.csv" % rprNextId, rprCands );
        # echo "mv candidates_it${rprNextId}_tmp.csv candidates_it${rprNextId}.csv"
        # mv candidates_it${rprNextId}_tmp.csv candidates_it${rprNextId}.csv
        if os.path.isfile( "candidates_it%d_tmp.csv" % rprNextId ):
            shutil.move( "candidates_it%d_tmp.csv" % rprNextId, "candidates_it%d.csv" % rprNextId ); # move back tmp

    # Show candidates
    #my_exec "../globOptVis --show$flag3D --scale $scale --pop-limit $poplimit -p $rprCands -a $rprReprAssoc --title \"GlobOpt-repr_candidates\" $visdefparam &"

    # Formulate
    # my_exec "$executable --formulate$flag3D $formParams --scale $scale --cloud cloud.ply --unary $unary --pw $rprPw --cmp $cmp --constr-mode patch --dir-bias $dirbias --patch-pop-limit $poplimit --angle-gens $anglegens --candidates $rprCands -a $rprReprAssoc --freq-weight $freqweight --cost-fn $pwCostFunc"
    cmd = "%s --formulate%s --scale %f --unary %f --pw %f --spat-weight %f --spat-dist-mult 2. --patch-pop-limit %d --angle-gens %s --cloud %s --candidates %s -a %s --collapse-angle-deg %f --trunc-angle %f --constr-mode patch --dir-bias 0 --no-clusters --cmp 0 --freq-weight 0 --cost-fn spatsqrt" \
           % ( args.rapterExec, args.flag3D, args.scale, args.data, rprPw, args.spatial, args.popLimit, angleGens, args.cloud, rprCands, rprReprAssoc, collapseThreshDeg, args.angleLimit)
    call( cmd, args.dry );

    rprDiagF    = "diag_it%d.gv" % rprIter;
    rprDiagFTmp = "%s%s" % (rprDiagF,"RprTmp");
    if not args.dry:
        # echo "cp primitives_it${rprIter}.bonmin.csv primitives_it${rprIter}_rprtmp.csv"
        # cp primitives_it${rprIter}.bonmin.csv primitives_it${rprIter}_rprtmp.csv
        shutil.copyfile( "primitives_it%d.bonmin.csv" % rprIter, "primitives_it%d_rprtmp.csv" % rprIter );
        if os.path.isfile( rprDiagF ): # backup diag_itx.gv
            #echo "mv $rprDiagF $rprDiagFTmp";
            # mv $rprDiagF "$rprDiagFTmp"
            shutil.move( rprDiagF, rprDiagFTmp );

    # my_exec "$executable --solver$flag3D bonmin --problem problem -v --time -1 --angle-gens $anglegens --bmode $algCode --candidates $rprCands"
    cmd = "%s --solver%s bonmin --problem problem -v --time -1 --angle-gens %s --bmode %d --candidates %s" \
        % (args.rapterExec, args.flag3D, angleGens, args.algCode, rprCands )
    call (cmd, args.dry );

    if not args.dry:
        # echo "cp primitives_it${rprIter}.bonmin.csv $rprReprOpt"
        # cp primitives_it${rprIter}.bonmin.csv $rprReprOpt
        shutil.copyfile( "primitives_it%d.bonmin.csv" % rprIter, rprReprOpt );

        # echo "cp primitives_it${rprIter}_rprtmp.csv primitives_it${rprIter}.bonmin.csv"
        # cp primitives_it${rprIter}_rprtmp.csv primitives_it${rprIter}.bonmin.csv
        shutil.copyfile( "primitives_it%d_rprtmp.csv" % rprIter, "primitives_it%d.bonmin.csv" % rprIter );

        # echo "mv $rprDiagF diag_it${rprIter}.lvl2.gv"
        # mv $rprDiagF diag_it${rprIter}.lvl2.gv
        shutil.move( rprDiagF, "diag_it%d.lvl2.gv" % rprIter );

        # restore diag_itx.gv
        if os.path.isfile( rprDiagFTmp ):
            # echo "mv $rprDiagFTmp $rprDiagF"
            # mv "$rprDiagFTmp" $rprDiagF
            shutil.move( rprDiagFTmp, rprDiagF );
            # rm "$rprDiagFTmp";
            #os.remove( rprDiagFTmp );

    #my_exec "../globOptVis --show$flag3D -p $rprReprOpt -a $rprReprAssoc --title \"GlobOpt-RepresentativesOptimized\" --scale $scale --pop-limit $poplimit $visdefparam &"

    # apply representatives - outputs subs.csv
    #my_exec "$executable --representBack$flag3D --repr $rprReprOpt -p $rprPrims -a $rprAssoc -sc $scale --cloud cloud.ply --angle-gens $anglegens"
    cmd = "%s --representBack%s --repr %s -p %s -a %s -sc %f --cloud %s --angle-gens %s" \
           % (args.rapterExec, args.flag3D, rprReprOpt, rprPrims, rprAssoc, args.scale, args.cloud, angleGens );
    call( cmd, args.dry );

    if not args.dry:
        # echo "mv $rprPrims $rprPrimBak"
        # mv $rprPrims $rprPrimBak
        shutil.move( rprPrims, rprPrimBak );
        # echo "mv subs.csv $rprPrims" #substitue for input
        # mv subs.csv $rprPrims
        shutil.move( "subs.csv", rprPrims );


parser = argparse.ArgumentParser()

suggestedGroup = parser.add_argument_group('suggested');
suggestedGroup.add_argument( "-s"  , "--scale"      , dest="scale"        , type=float, default=0.05, help="Scale (rho) parameter, the smallest feature size to preserve [0.001..0.05]", required=True)
suggestedGroup.add_argument( "--al", "--angle-limit", dest="angleLimit"   , type=float, default=15  , help="Angle threshlod (tau) parameter in degrees [5..45]")
suggestedGroup.add_argument( "--pw", "--pairwise"   , dest="pw"           , type=float, default=1.0 , help="Weight of pairwise term [0.1..10^6]" )
suggestedGroup.add_argument( "-t" , "--area-thresh-start"  , dest="smallThreshMult", type=float, default= 4., help="Start with planes, that are scale * smallThreshMult large. Increase this, if optimisation too slow. [powers of 2].")

optionalGroup = parser.add_argument_group('optional');
optionalGroup.add_argument( "--ag", "--angle-gens" , dest="angleGens"    , type=float, default=[0,90], help="Weight of pairwise term [0.1..10^6]", action="append" )
optionalGroup.add_argument( "--it", "--iterations" , dest="nbExtraIterations", type=int, default=15, help="How many iterations to run [5..20]")
optionalGroup.add_argument( "--cl", "--cloud"      , dest="cloud"        , type=str  , default = "cloud.ply", help="Pointcloud in ply format [cloud.ply]");
optionalGroup.add_argument( "-l"  , "--lines"      , dest="lines"        , action="store_true"     , help="Work in 2D with lines instead of planes." )

runOptGroup = parser.add_argument_group('run options');
runOptGroup.add_argument( "--dry", action="store_true"                 , help="Show the calls, but don't run." )
runOptGroup.add_argument( "--no-vis", dest="noVis", action="store_false", default = False, help="Disable visualization (enabled by default)" )

optionalGroup.add_argument( "--pl", "--popLimit"   , dest="popLimit"     , type=int  , default=5   , help="Filters primitives having less than this many points assigned [3..100]")
optionalGroup.add_argument( "--sp", "--spatial"    , dest="spatial"      , type=float,               help="Weight of spatial term [0.1, pw/10., pw/5., pw/2.]" )
optionalGroup.add_argument("--vl" , "--var-limit"  , dest="variableLimit", type=int  , default=1000, help="Maximum number of variables (primitives) for the optimisation. [500..3000]")

optionalGroup.add_argument( "-d"  , "--data"       , dest="data"         , type=float, default=1e5 , help="Weight of data term [10^5, 10^6]" )
optionalGroup.add_argument( "-p"  , "--primitives" , dest="primitives"   , type=str  ,               help="Input primitives, e.g. existing segmentation segments.csv" )
optionalGroup.add_argument( "-a"  , "--assoc"      , dest="associations" , type=str  ,               help="Input point-primitive associations, e.g. existing segmentation's points_segments.csv" )


optionalGroup.add_argument("--segment-scale-mult"          , dest="segmentScaleMultiplier", type=float, default=1.0, help="Multiply scale by this value for the segmentation step. [0.5, 1.0, 2.0]")
optionalGroup.add_argument("--ald", "--angle-limit-divisor", dest="angleLimitDivisor"     , type=float, default=1.0, help="Divide angle threshold (tau) by this number for candidate generation. [2.0, 1.0, 0.5]")
optionalGroup.add_argument("--alg-code"                    , dest="algCode"               , type=int  , default=0  , help="Bonmin algorithm enum codes. 0: B_BB, 1: OA, 2: QG, 3: Hyb, 4: ECP, 5: IFP. [0]");

args = parser.parse_args()

if not os.path.isfile(args.cloud):
    print("Need \"%s\" to exist, assuming it's the pointcloud" % args.cloud );
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

setattr( args, "rapterExec", rapterExec );
if not args.lines:
    setattr( args, "flag3D"     ,"3D"             );
    setattr( args, "tripletSafe","--triplet-safe" ); # Values: ["", "--triplet-safe"]
    useAllGens = min(5, args.nbExtraIterations-1 );                    # start with parallel generation only
else:
    setattr( args, "flag3D"     ,"" );
    setattr( args, "tripletSafe","" );
    useAllGens = 0

########################################################################################################################
# Do segmentation
if not args.primitives or not args.associaitons:
    cmd = "%s --segment%s --scale %f --angle-limit %f --angle-gens %s --patch-pop-limit %d --dist-limit-mult %f --cloud %s" \
          % ( rapterExec, args.flag3D, args.scale, args.angleLimit, args.angleGensStr, args.popLimit, args.segmentScaleMultiplier, args.cloud )
    call( cmd, args.dry );
    # # save output
    if ( os.path.isfile("patches.csv") and os.path.isfile("points_primitives.csv") ):
        if os.path.isfile("segments.csv"):
            shutil.copyfile( "segments.csv", "segments.csv.bak" );
        shutil.copyfile( "patches.csv", "segments.csv" )
        if os.path.isfile("points_segments.csv"):
            shutil.copyfile( "points_segments.csv", "points_segments.csv.bak" );
        shutil.copyfile( "points_primitives.csv", "points_segments.csv" )

    if not args.noVis:
        show( "segments.csv", "points_segments.csv", "\"RAPter - Segmentation\"", args );

########################################################################################################################

########################################################################################################################

angleGens           = "0"
candAngleGens       = "0"                       # used to mirror anglegens, but keep const "0" for generate
primitives          = "patches.csv"
associations        = "points_primitives.csv"
keepSingles         = "--keep-singles"
allowPromoted       = "--allow-promoted"
smallThreshDiv      = 2.                        # area threshold stepsize
smallThreshLimit    = 0.                        # when to stop decreasing area threshold
promRem             = 0                         # remaining primitives to promote
collapseThreshDeg   = 0.4                       # initialize optimisation with the closest two orientations merged, if their difference is < colleapseThreshDeg degrees.
adopt               = 0
adoptChanged        = False
decreaseLevel       = False

iteration = 0

while iteration <= args.nbExtraIterations:

    # decresase, unless there is more to do on the same level
    if decreaseLevel:
        args.smallThreshMult = float(int(args.smallThreshMult / smallThreshDiv));

    # if we reached the bottom working scale (ideally 0)
    if args.smallThreshMult <= smallThreshLimit:
        args.smallThreshMult = int(smallThreshLimit) # make sure it's integer
        if decreaseLevel:                       # if we don't have to promote any more patches on this level
            adopt = "1"                         # if we promoted all patches, we can allow points to get re-assigned
            if not adoptChanged:
                adoptChanged = True             # only enter here once
                useAllGens   = iteration + 2    # if we promoted all patches in the scene, do a 90 round
                args.nbExtraIterations = max(args.nbExtraIterations,useAllGens + 3)   # do k more rounds after the 90 round

    # reset to false, meaning we will continue decreasing, unless generate flips it again
    decreaseLevel = True

    print( "smallThreshMult: %d" % args.smallThreshMult );
    print( "__________________________________________________________" );
    print( "Start iteration %d" % iteration );

    prevId = iteration - 1;
    nextId = iteration + 1;

    if iteration > 0:
        primitives   = "primitives_merged_it%d.csv" % prevId;
        associations = "points_primitives_it%d.csv" % prevId;

    # (2) Generate - generate candidates
    cmd = "%s --generate%s -sc %f -al %f -ald %f --patch-pop-limit %d -p %s --assoc %s --cloud %s --angle-gens %s --small-thresh-mult %f --var-limit %d --small-mode 0 %s %s %s" \
            % (rapterExec, args.flag3D, args.scale, args.angleLimit, args.angleLimitDivisor, args.popLimit, primitives, associations, args.cloud, candAngleGens, \
               args.smallThreshMult, args.variableLimit, \
               args.tripletSafe, keepSingles, allowPromoted );
    promRem = call( cmd, args.dry, True );
    print "[rapter.py] Remaining smalls to promote: ", promRem
    # Don't decrase area threshold until no more candidates to promote
    if promRem != 0:
        decreaseLevel = False; # set to true each iteration

    # (3) Formulate - create optimisation problem
    cmd = "%s --formulate%s --scale %f --unary %f --pw %f --spat-weight %f --spat-dist-mult 2. --patch-pop-limit %d --angle-gens %s --cloud %s --candidates candidates_it%d.csv -a %s --collapse-angle-deg %f --trunc-angle %f --constr-mode patch --dir-bias 0 --no-clusters --cmp 0 --freq-weight 0 --cost-fn spatsqrt" \
           % ( rapterExec, args.flag3D, args.scale, args.data, args.pw, args.spatial, args.popLimit, angleGens, args.cloud, iteration, associations, collapseThreshDeg, args.angleLimit)
    call( cmd, args.dry );

    # (4) Solve
    cmd = "%s --solver%s bonmin --problem problem -v --time -1 --bmode %d --angle-gens %s --candidates candidates_it%d.csv --cloud %s" \
           % ( rapterExec, args.flag3D, args.algCode, angleGens, iteration, args.cloud );
    call( cmd, args.dry )

    if not args.noVis:
        show( "primitives_it%d.bonmin.csv" % iteration, associations, "\"RAPter - Iteration%d\"" % iteration, args );

    if iteration == useAllGens:
        angleGens     = ','.join( str(e) for e in args.angleGens );
        candAngleGens = angleGens;

    # TODO
    runRepr( "primitives_it%d.bonmin.csv" % iteration, associations, iteration, args, angleGens, keepSingles );

    # (6) CoPlanarity
    cmd = "%s --merge%s --scale %f --adopt %s --prims primitives_it%d.bonmin.csv -a %s --angle-gens %s --patch-pop-limit %d --cloud %s" \
          % ( rapterExec, args.flag3D, args.scale, adopt, iteration, associations, angleGens, args.popLimit, args.cloud );
    call( cmd, args.dry )

    # Don't copy promoted patches' directions to other patches after 4 iterations (c==3), since they are not reliable anymore
    if iteration == 3:
        allowPromoted = ""

    # Don't throw away single directions before the 3rd (c==1) iteration.
    # This will keep large patches, even if they don't copy to anywhere for later.
    if iteration == 1:
        keepSingles = ""

    # If we are still promoting small patches on this working scale, make sure to run more iterations
    if iteration == args.nbExtraIterations and promRem != 0:
        args.nbExtraIterations += 1;

    # Increment iteration counter (while loop)
    iteration += 1;
