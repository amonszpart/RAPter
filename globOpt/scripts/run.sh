# Executable path
executable="../glob_opt";
execLog="lastRun.log"
execPyRelGraph="python ../readGraphProperties.py"
execPyStats="python ../collectStatistics.py"
###############################################
doExecPy=true           # Set to true if generating graphs
dryRun=false            # Values: [true/false]. If "true", don't execute, just show calls
noVis=false             # Only show premerge and final output

# startAt=0: do segment
# startAt=1: don't do segment, do premerge
# startAt=2: don't do premerge, do iteration 0
startAt=0
nbExtraIter=40          # Values: [ 1..20 ] iteration count. Default: 2.
reprPwMult=3            # Values: [ 1, 5 ,10]. Multiplies pw with this number when doing lvl2 (representatives)
use90=24                # Values: [ 100 ] Will be re-set by runscript
extendedAngleGens="90"  # Values: [ "0", "90", ... ]

anglegens="0"           # Values: ["0", "60", "90", "60,90", ... ] Desired angle generators in degrees. Default: 90.
unary=100000            # Values: [1000, 1000000]
spatWeight=0            # Values: [ 10, 0 ] (Penalty that's added, if two primitives with different directions are closer than 2 x scale)
truncAngle=$anglelimit  # Values: [ 0, $anglelimit, 0.15 ] (Pairwise cost truncation angle in radians)
smallThreshDiv="2"      # Values: [ 2, 3, ... ] Stepsize of working scale.

# Parse runScripts
runFileName=`basename $0`;
runPath=$(dirname "`file $0 | awk '{print $5}' |cut -c 2-`");
echo "runPath:" $runPath
source "$runPath/runRepr.sh"; echo "loaded runRepr.sh"
source "$runPath/runFuncs.sh"; echo "loaded runFuncs.sh"

function print_usage() {
        echo "usage:\t run.sh scale anglelimit pairwisecost [pop-limit] [3D] [smallThreshStart]"
        echo "example:\t run.sh 0.03 0.2 1 25 3D 64"
        echo "showPearl: ../globOptVis --show --scale 0.05 --pop-limit 0 -p primitives.pearl.csv -a points_primitives.pearl.csv --title \"Pearl\" --use-tags --no-clusters --no-pop"
}

# parse scale
if [[ -z "$1" ]]; then print_usage; exit 1; else scale=$1; fi
# parse angle-limit
if [[ -z "$2" ]]; then print_usage; exit 1; else anglelimit=$2; fi
# parse pairwise cost
if [[ -z "$3" ]]; then print_usage; exit 1; else pw=$3; fi
# parse population limit
if [ -n "$4" ]; then poplimit=$4; else poplimit=5; fi
# parse 3D
if [ -n "$5" ]; then flag3D=$5; else flag3D=""; fi
# parse smallThresh (a size(spatialSignificance) threshold for primitives to work with)
if [ -n "$6" ]; then smallThresh=$6; else smallThresh="1"; fi

########################################
#### We usually don't change the following:
# Moved to argument 6
#smallThresh="256" # smallThresh * scale is the small threshold # 5 was good for most of the stuff, except big scenes (kinect, lanslevillard)
smallThreshlimit="0"    # last threshold

dirbias="0";	# Values: [ 0, *1* ] (Not-same-dir-id cost offset. Default: 0. Don't use, if freqweight is on.)
freqweight="0"; # Values: [ 0, *1000* ] (Dataterm = (freqweight / #instances) * datacost)
adopt="0";      # Values: [ *0*, 1] (Adopt points argument. Default: 0)
cmp=0           # Values: [0,1]
allowPromoted="--allow-promoted" # Will be turned off after 3 iterations at the end of loop
keepSingles="--keep-singles"     # Will be turned off after 1 iteration at the end of loop

noClusters="--no-clusters"  # Values: [ "", "--no-clusters" ] (Flag that turns spatial clustering extra variables off)
useAngleGen=""              # Values: [ "", "--use-angle-gen" ] (Flag which makes not same directioned primitives not to penalize eachother)
formParams="$noClusters $useAngleGen --spat-weight $spatWeight --trunc-angle $truncAngle"

# In candidate generation, divide angle limit with this to match copies. Default: 1. Set to 10, if too many candidates (variables).
cand_anglediv="1";          # for 3D: "2.5";
# Multiply scale by this number to get the segmentation (regionGrowing) spatial distance
segmentScaleMultiplier="1"; # for 3D: "2.5";
pwCostFunc="spatsqrt"       # Spatial cost function.

visdefparam="--angle-gens $anglegens --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-rel --no-scale --bg-colour .9,.9,.9" #"--use-tags --no-clusters" #--ids
iterationConstrMode="patch" # what to add in the second iteration formulate. Default: 0 (everyPatchNeedsDirection), experimental: 2 (largePatchesNeedDirection).

safeMode="";                #"--safe-mode"; #"--safe-mode" # "--safe-mode" for new, or "" for old version
variableLimit=1100; # 1300; # Safe mode gets turned on, and generate rerun, if candidates exceed this number (1300)
premerge=1                  # Call merge after segmentation 0/1
algCode=0                   # 0==B_BB, OA, QG, 3==Hyb, ECP, IFP
#######################################
#######################################

echo "scale: $scale"
echo "pw: $pw"
echo "pop-limit: $poplimit"
echo "anglegens: $anglegens"
echo "angle-limit: $anglelimit"
echo "freq-weight: $freqweight"
echo "region-grow-scale-mult: $segmentScaleMultiplier"
echo "smallThresh: " $smallThresh
echo "pw cost: " $pwCostFunc

rm $execLog

# call it
save_args $0 $@ "--reprPwMult $reprPwMult --extendedAngleGens $extendedAngleGens $formParams --freqweight" $freqweight "--angle-limit" $anglelimit "--segment-scale-mult" $segmentScaleMultiplier "--adopt" $adopt "--dirbias" $dirbias "--cand-anglediv" ${cand_anglediv} "--angle-gens" $anglegens "--cost-fn" $pwCostFunc "--small-thresh" $smallThresh "--small-thresh-limit" $smallThreshLimit "--smallThreshDiv" $smallThreshDiv --unary $unary --cmp $cmp

#####
# Check if the gt folder exists, in that case compute the primitive comparisons
correspondance=false
correspondance_exe="../corresp"
correspondance_gtprim="./gt/primitives.csv"
correspondance_gtassing="./gt/points_primitives.csv"
if [ -e "$correspondance_gtprim" ]; then
  echo "Ground truth folder detected, compute primitive correspondances" 
  correspondance=true  
fi
####

# Backup energy.csv
mv energy.csv energy.csv.bak

input="patches.csv";
assoc="points_primitives.csv";

# [0] Segmentation. OUTPUT: patches.csv, points_primitives.csv
if [ $startAt -eq 0 ]; then
    my_exec "$executable --segment$flag3D --angle-limit $anglelimit --scale $scale --dist-limit-mult $segmentScaleMultiplier --angle-gens $anglegens"
    # save output
    cp $input "segments.csv"
    cp $assoc "points_segments.csv"
fi

# PreMerge?
if [ $startAt -le 1 ]; then
    # merge before start
    if [ $premerge -ne 0 ]; then
        echo "performing premerge"
        my_exec "$executable --merge$flag3D --scale $scale --adopt $adopt --prims $input -a $assoc --angle-gens $anglegens --patch-pop-limit $poplimit"

        # save patches
        cp patches.csv_merged_it-1.csv $input
        cp points_primitives_it-1.csv $assoc
        my_exec "../globOptVis --show$flag3D --scale $scale --use-tags --pop-limit $poplimit -p patches.csv -a $assoc --normals 100 --title \"GlobOpt - PreMerge output\" --no-clusters --no-pop --no-rel --bg-colour .9,.9,.9 --angle-gens $anglegens &"
        if [ "$noVis" = false ] ; then
            my_exec "../globOptVis --show$flag3D --scale $scale --pop-limit $poplimit -p segments.csv -a points_segments.csv --title \"GlobOpt - Segmentation output\" $visdefparam --dir-colours --no-rel &"
        fi
    fi #...if premerge
fi #...startAt <= 1

# If true, stay on the same level for one more iteration (too many variables)
decrease_level=false;
adoptChanged=false;

c=0         # Iteration count variable, goes from 0 .. nbExtraIter
promRem=0   # Count of remaining patches to promote
while [ $c -le $nbExtraIter ];
do
    # decresase, unless there is more to do on the same level
    if $decrease_level; then
        smallThresh=`../divide.py $smallThresh $smallThreshDiv`;
        smallThresh=$smallThresh;
    fi
    
    # if we reached the bottom working scale (ideally 0)
    if [ $smallThresh -lt $smallThreshlimit ] || [ $smallThresh -eq "0" ]; then
        smallThresh=$smallThreshlimit           # make sure it's integer
        if $decrease_level; then                # if we don't have to promote any more patches on this level
            adopt="1"                           # if we promoted all patches, we can allow points to get re-assigned
            if ! $adoptChanged; then
                use90=$(( $c + 1 ));            # if we promoted all patches in the scene, do a 90 round
                nbExtraIter=$(( $use90 + 3 ))   # do k more rounds after the 90 round
                adoptChanged=true               # only enter here once
                # debug:
                echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<< use90: $use90, nbExtraIter: $nbExtraIter, c: $c, \
                      smallThresh: $smallThresh, decrl: $decrease_level, promRem: $promRem"
            fi
        fi
    fi

    # reset to false, meaning we will continue decreasing, unless generate flips it again
    decrease_level=true

    echo "smallThreshMult: " $smallThresh
	  echo "__________________________________________________________";
    echo "Start iteration $c"
    
    prevId=`expr $c - 1`;
    nextId=`expr $c + 1`;

    if [ $c -gt 0 ]; then
        input="primitives_merged_it$prevId.csv";
        assoc="points_primitives_it$prevId.csv";
    fi

    if [ $c -ge $(($startAt - 2)) ]; then
        # Generate candidates from output of first. OUT: candidates_it$c.csv. #small-mode 2: small patches receive all candidates
        my_exec2 "$executable --generate$flag3D $allowPromoted $keepSingles -sc $scale -al $anglelimit -ald ${cand_anglediv} --small-mode 0 --patch-pop-limit $poplimit -p $input --assoc $assoc --angle-gens $anglegens --small-thresh-mult $smallThresh --var-limit $variableLimit $safeMode"
        echo "Remaining smalls to promote: $myresult"
        # count of remaining patches to promote
        promRem=$myresult 

        # don't decrease spatial threshold, if there are still variables to promote
        if [ $myresult -ne "0" ]; then decrease_level=false; fi

        # Formulate optimization problem. OUT: "problem" directory. --constr-mode 2: largePatchesNeedDirectionConstraint
        my_exec "$executable --formulate$flag3D --scale $scale --cloud cloud.ply --unary $unary --pw $pw --cmp $cmp --constr-mode $iterationConstrMode --dir-bias $dirbias --patch-pop-limit $poplimit --angle-gens $anglegens --candidates candidates_it$c.csv -a $assoc --freq-weight $freqweight  --cost-fn $pwCostFunc $formParams"

        # Solve optimization problem. OUT: primitives_it$c.bonmin.csv
        my_exec "$executable --solver$flag3D bonmin --problem problem -v --time -1 --bmode $algCode --angle-gens $anglegens --candidates candidates_it$c.csv"
        
        if [ "$correspondance" = true ] ; then
            my_exec "$correspondance_exe  $correspondance_gtprim $correspondance_gtassing primitives_it$c.bonmin.csv $assoc cloud.ply $scale"
        fi

        # Show output
        if [ "$noVis" = false ] ; then
            my_exec "../globOptVis --show$flag3D --scale $scale --pop-limit $poplimit -p primitives_it$c.bonmin.csv -a $assoc --title \"GlobOpt - $c iteration output\" $visdefparam &"
        fi
        
        # Generate relation graphs
        if [ "$doExecPy" = true ]; then
            my_exec "$execPyRelGraph primitives_it$c.bonmin.csv $assoc cloud.ply --angles $anglegens --iteration $c"
        fi

        # Use 90 degrees only for a single iteration's lvl2
        if [ $c -eq $use90 ]; then
            _bakAngleGens=$anglegens;
            anglegens=$extendedAngleGens;
        fi

        # Representatives (lvl2)
        echo "[call] runRepr primitives_it$c.bonmin.csv $assoc $c"
        runRepr primitives_it$c.bonmin.csv $assoc $c

        #if [ $(( $c - 1 )) -eq $use90 ]; then
        if [ $(( $c )) -eq $use90 ]; then
            anglegens=$_bakAngleGens;
        fi

        # Merge/show parallel
        if [ $c -lt $nbExtraIter ]; then
            # Merge adjacent candidates with same dir id. OUT: primitives_merged_it$c.csv, points_primitives_it$c.csv
            my_exec "$executable --merge$flag3D --scale $scale --adopt $adopt --prims primitives_it$c.bonmin.csv -a $assoc --angle-gens $anglegens --patch-pop-limit $poplimit"
        else
            my_exec "../globOptVis --show$flag3D --scale $scale --pop-limit $poplimit -p primitives_it$c.bonmin.csv -a $assoc --title \"GlobOpt - [Dir-Colours] $c iteration output\" $visdefparam --paral-colours --no-rel &"
        fi

        #my_exec "$executable --energy --formulate$flag3D --scale $scale --cloud cloud.ply --unary $unary --pw $pw --cmp $cmp --constr-mode $iterationConstrMode --dir-bias $dirbias --patch-pop-limit $poplimit --angle-gens $anglegens --candidates primitives_it$c.bonmin.csv -a $assoc --freq-weight $freqweight  --cost-fn $pwCostFunc $formParams"
    fi

    # Don't copy promoted patches' directions to other patches after 4 iterations (c==3), since they are not reliable anymore
    if [ $c -eq 3 ]; then allowPromoted=""; fi
    # Don't throw away single directions before the 3rd (c==1) iteration. 
    # This will keep large patches, even if they don't copy to anywhere for later.
    if [ $c -eq 1 ]; then keepSingles=""; fi

    # If we are still promoting small patches on this working scale, make sure to run more iterations
    if [ $c -eq $nbExtraIter ] && [ $promRem -ne "0" ]; then
        nbExtraIter=$(( $nbExtraIter + 1 ))
    fi

    # Increment iteration counter (while loop)
    c=$(( $c + 1 ))
done

if [ "$doExecPy" = true ]; then
    my_exec "$execPyStats .  --angles $anglegens"
fi

#energies
# ../pearl --scale 0.05 --cloud cloud.ply --prims patches.csv --assoc points_primitives.csv --pw 1000 --cmp 1000

# remove safemode
# put back startat
# put back adopt?
