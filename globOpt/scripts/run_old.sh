#    Executable path
executable="../glob_opt";
execLog="lastRun.log"
execPyRelGraph="python ../readGraphProperties.py"
execPyStats="python ../collectStatistics.py"
doExecPy=false # set to true if generating graphs

###############################################

function print_usage() {
        echo "usage:\t run.sh scale anglelimit pairwisecost [pop-limit] [3D] [smallThreshStart]"
        echo "example:\t run.sh 0.03 0.2 1 25 3D 64"
        echo "showPearl: ../globOptVis --show --scale 0.05 --pop-limit 0 -p primitives.pearl.csv -a points_primitives.pearl.csv --title \"Pearl\" --use-tags --no-clusters --no-pop"
}

# parse scale
if [[ -z "$1" ]]; then
	print_usage
	exit 1
else
	scale=$1;
fi

# parse angle-limit
if [[ -z "$2" ]]; then
	print_usage
	exit 1
else
	anglelimit=$2
fi

# parse pairwise cost
if [[ -z "$3" ]]; then
	print_usage
	exit 1
else
	pw=$3
fi

# parse population limit
if [ -n "$4" ]; then
	poplimit=$4
else
	poplimit=5
fi

# parse 3D
if [ -n "$5" ]; then
	flag3D=$5;
else
	flag3D="";
fi

#parse smallThresh (a size(spatialSignificance) threshold for primitives to work with)
if [ -n "$6" ]; then
        smallThresh=$6;
else
        smallThresh="1";
fi

anglegens="0"; # Values: ["60", "90", "60,90", ... ] Desired angle generators in degrees. Default: 90.
nbExtraIter=14;  # Values: [ 1..20 ] iteration count. Default: 2.
dirbias="0";	# Values: [ 0, *1* ] (Not-same-dir-id cost offset. Default: 0. Don't use, if freqweight is on.)
freqweight="0"; # Values: [ 0, *1000* ] (Dataterm = (freqweight / #instances) * datacost)
adopt="0";      # Values: [ *0*, 1] (Adopt points argument. Default: 0)
# In candidate generation, divide angle limit with this to match copies. Default: 1. Set to 10, if too many candidates (variables).
cand_anglediv="1";# for 3D: "2.5";
# multiply scale by this number to get the segmentation (regionGrowing) spatial distance
segmentScaleMultiplier="1";# for 3D: "2.5";
pwCostFunc="spatsqrt" # spatial cost function. TODO: reactivate sqrt (does not compile for now)
unary=1000 #1000000

noClusters="--no-clusters" # Values: [ "", "--no-clusters" ] (Flag that turns spatial clustering extra variables off)
useAngleGen="" # Values: [ "", "--use-angle-gen" ] (Flag which makes not same directioned primitives not to penalize eachother)
spatWeight=0 # Values: [ 10, 0 ] (Penalty that's added, if two primitives with different directions are closer than 2 x scale)
truncAngle=$anglelimit # Values: [ 0, $anglelimit, 0.15 ] (Pairwise cost truncation angle in radians)
formParams="$noClusters $useAngleGen --spat-weight $spatWeight --trunc-angle $truncAngle"

visdefparam="--angle-gens $anglegens --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-rel --no-scale --bg-colour .9,.9,.9" #"--use-tags --no-clusters" #--ids
firstConstrMode="patch" # what to add in the first run formulate. Default: 0 (everyPatchNeedsDirection), experimental: 2 (largePatchesNeedDirection).
iterationConstrMode="patch" # what to add in the second iteration formulate. Default: 0 (everyPatchNeedsDirection), experimental: 2 (largePatchesNeedDirection).

# startAt=0: do segment
# startAt=1: don't do segment, do iteration 0
# startAt=2: don't do iteration 0, do iteration 1
startAt=0

#moved to argument 6
#smallThresh="256" # smallThresh * scale is the small threshold # 5 was good for most of the stuff, except big scenes (kinect, lanslevillard)
#smallThresh="128"
smallThreshlimit="0" #last threshold
smallThreshDiv="2"; #stepsize

safeMode=""; #"--safe-mode"; #"--safe-mode" # "--safe-mode" for new, or "" for old version
variableLimit=1100; # 1300; # Safe mode gets turned on, and generate rerun, if candidates exceed this number (1300)
premerge=1 # call merge after segmentation 0/1
algCode=0 # 0==B_BB, OA, QG, 3==Hyb, ECP, IFP

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

# save run.log arguments to "run.log"
function save_args() {
	logfile="run.log"
	args=("$@") 
	# get number of elements 
	ELEMENTS=${#args[@]} 

	echo -n "[$(date +%D\ %T)] " >> $logfile
	# echo each element in array  
	# for loop 
	for (( i=0;i<$ELEMENTS;i++)); do 
	    echo -n "${args[${i}]} " >> $logfile
	done
	# endline
	echo -e -n "\n" >> $logfile
}

# call it
save_args $0 $@ "$formParams --freqweight" $freqweight "--angle-limit" $anglelimit "--segment-scale-mult" $segmentScaleMultiplier "--adopt" $adopt "--dirbias" $dirbias "--cand-anglediv" ${cand_anglediv} "--angle-gens" $anglegens "--cost-fn" $pwCostFunc "--small-thresh" $smallThresh "--small-thresh-limit" $smallThreshLimit "--smallThreshDiv" $smallThreshDiv

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

# show command before run
# stop script when the command fails
function my_exec() {
	echo "__________________________________________________________";
	echo -e "\n\n[CALLING] $1";
        echo -e $1"\n" >>$execLog
        eval $1;
  	if [ "$?" -ne "0" ]; then
	    echo "Error detected ($?). ABORT."
	    exit 1
  	fi
}

# won't halt, just report return value
function my_exec2() {
        echo "__________________________________________________________";
        echo -e "\n\n[CALLING] $1";
        echo -e $1"\n" >>$execLog
        eval $1;
        ret=$?;
        if [ "$ret" -ne "0" ]; then
            echo "Call returned $ret"
        fi

        myresult=$ret;
}

function countLines()
{
    echo `wc -l $1 | cut -f1 -d' '`;
}

mv energy.csv energy.csv.bak

# [0] Segmentation. OUTPUT: patches.csv, points_primitives.csv
if [ $startAt -eq 0 ]; then
    my_exec "$executable --segment$flag3D --angle-limit $anglelimit --scale $scale --dist-limit-mult $segmentScaleMultiplier --angle-gens $anglegens"
fi

input="patches.csv";
assoc="points_primitives.csv";

# show segment output
#my_exec "../globOptVis --show$flag3D --scale $scale --use-tags --pop-limit $poplimit -p patches.csv -a $assoc --normals 100  --no-clusters --title \"GlobOpt - Segment output\" --no-pop --no-rel &"

if [ $startAt -le 1 ]; then

# merge before start
if [ $premerge -ne 0 ]; then
    echo "performing premerge"
    my_exec "$executable --merge$flag3D --scale $scale --adopt $adopt --prims $input -a $assoc --angle-gens $anglegens --patch-pop-limit $poplimit"
    premergeLinesCnt=`wc -l patches.csv_merged_it-1.csv | cut -f1 -d' '`;
    patchesLinesCnt=`wc -l patches.csv | cut -f1 -d' '`;
    echo ">>>>>> PREMERGE produced  $premergeLinesCnt patches instead of $patchesLinesCnt";
    # save patches
    cp $input segments.csv
    cp $assoc points_segments.csv
    cp patches.csv_merged_it-1.csv $input
    cp points_primitives_it-1.csv $assoc
    my_exec "../globOptVis --show$flag3D --scale $scale --use-tags --pop-limit $poplimit -p patches.csv -a $assoc --normals 100 --title \"GlobOpt - PreMerge output\" --no-clusters --no-pop --no-rel --bg-colour .9,.9,.9 --angle-gens $anglegens &"
    my_exec "../globOptVis --show$flag3D --scale $scale --pop-limit $poplimit -p segments.csv -a points_segments.csv --title \"GlobOpt - Segmentation output\" $visdefparam --dir-colours --no-rel &"
fi

# Generate candidates. OUT: candidates_it0.csv. #small-mode : small patches don't receive any candidates
my_exec2 "$executable --generate$flag3D -sc $scale -al $anglelimit -ald ${cand_anglediv} --small-mode 0 --patch-pop-limit $poplimit -p $input --assoc $assoc --angle-gens $anglegens --small-thresh-mult $smallThresh"

# Formulate optimization problem. OUT: "problem" directory. constr-mode 2: largePatchesNeedDirectionConstraint
my_exec "$executable --formulate$flag3D --scale $scale --cloud cloud.ply --unary $unary --pw $pw --cmp 1 --constr-mode $firstConstrMode     --dir-bias $dirbias --patch-pop-limit $poplimit --angle-gens $anglegens --candidates candidates_it0.csv -a $assoc --freq-weight $freqweight --cost-fn $pwCostFunc $formParams"

# Solve optimization problem. OUT: primitives_it0.bonmin.csv
my_exec "$executable --solver$flag3D bonmin --problem problem -v --time -1 --bmode $algCode --angle-gens $anglegens --candidates candidates_it0.csv"

if [ "$correspondance" = true ] ; then
    my_exec "$correspondance_exe  $correspondance_gtprim $correspondance_gtassing primitives_it0.bonmin.csv $assoc cloud.ply $scale"
fi

# Show output of first iteration.
my_exec "../globOptVis --show$flag3D --scale $scale --pop-limit $poplimit -p primitives_it0.bonmin.csv  -a $assoc --title \"GlobOpt - 0th iteration output\" $visdefparam &"
# !! set flag3D to "3D" and recomment these two lines to work with 3D
#energies

# Merge adjacent candidates with same dir id. OUT: primitives_merged_it0.csv, points_primitives_it0.csv
my_exec "$executable --merge$flag3D --scale $scale --adopt $adopt --prims primitives_it0.bonmin.csv -a $assoc --angle-gens $anglegens --patch-pop-limit $poplimit"
#my_exec "$executable --merge$flag3D --scale $scale --adopt $adopt --prims primitives_it$c.bonmin.csv -a $assoc --angle-gens $anglegens --patch-pop-limit $poplimit"
#cp primitives_it0.bonmin.csv primitives_merged_it0.csv
#cp points_primitives.csv points_primitives_it0.csv

# Show output of first merge.
#my_exec "../globOptVis --show$flag3D --scale $scale --pop-limit $poplimit -p primitives_merged_it0.csv -a points_primitives_it0.csv --title \"GlobOpt - Merged 1st iteration output\" $visdefparam  &"
#my_exec "$executable --formulate$flag3D --energy --scale $scale --cloud cloud.ply --unary $unary --pw $pw --cmp 1 --constr-mode $firstConstrMode --dir-bias $dirbias --patch-pop-limit $poplimit --angle-gens $anglegens --candidates primitives_it0.bonmin.csv -a $assoc --freq-weight $freqweight --cost-fn $pwCostFunc $formParams"

# generate relation graphs
if [ "$doExecPy" = true ]; then
    my_exec "$execPyRelGraph primitives_it0.bonmin.csv points_primitives_it0.csv cloud.ply --angles $anglegens --iteration 0"
fi

fi #startAt <= 1

# if true, stay on the same level for one more iteration (too many variables)
decrease_level=true;
# converged is false, until two outputs are not distinguishable by diff
converged=false;

for c in $(seq 1 $nbExtraIter)
do
    # decresase, unless there is more to do on the same level
    if $decrease_level; then
        smallThresh=`../divide.py $smallThresh $smallThreshDiv`;
        smallThresh=$smallThresh;
        echo "!!!!!!!!!!! decreased !!!!!!!!!!!!";
    fi
    # reset to false, meaning we will continue decreasing, unless generate flips it again
    decrease_level=true
    
    if [ $smallThresh -lt $smallThreshlimit ] || [ $smallThresh -eq "0" ]; then
        smallThresh=$smallThreshlimit
        #if [ $startAt -eq 0 ]; then
            adopt="1"
        #fi
    fi

    echo "smallThreshMult: " $smallThresh
	  echo "__________________________________________________________";
    echo "Start iteration $c"
    
    prevId=`expr $c - 1`;
    nextId=`expr $c + 1`;

    input="primitives_merged_it$prevId.csv";
    assoc="points_primitives_it$prevId.csv";


    if [ $c -ge $(($startAt - 1)) ]; then
        # Generate candidates from output of first. OUT: candidates_it$c.csv. #small-mode 2: small patches receive all candidates
        my_exec2 "$executable --generate$flag3D -sc $scale -al $anglelimit -ald ${cand_anglediv} --small-mode 0 --patch-pop-limit $poplimit -p $input --assoc $assoc --angle-gens $anglegens --small-thresh-mult $smallThresh --var-limit $variableLimit $safeMode"
        echo "CNT: $myresult"
        if [ $myresult -ne "0" ]; then
            decrease_level=false;
        fi

        if $decrease_level; then
            echo "decrease_level: TRUE";
        else
            echo "decrease_level: FALSE";
        fi

        # Formulate optimization problem. OUT: "problem" directory. --constr-mode 2: largePatchesNeedDirectionConstraint
        my_exec "$executable --formulate$flag3D --scale $scale --cloud cloud.ply --unary $unary --pw $pw --cmp 1 --constr-mode $iterationConstrMode --dir-bias $dirbias --patch-pop-limit $poplimit --angle-gens $anglegens --candidates candidates_it$c.csv -a $assoc --freq-weight $freqweight  --cost-fn $pwCostFunc $formParams"

        # Solve optimization problem. OUT: primitives_it$c.bonmin.csv
        my_exec "$executable --solver$flag3D bonmin --problem problem -v --time -1 --bmode $algCode --angle-gens $anglegens --candidates candidates_it$c.csv"
        
        if [ "$correspondance" = true ] ; then
            my_exec "$correspondance_exe  $correspondance_gtprim $correspondance_gtassing primitives_it$c.bonmin.csv $assoc cloud.ply $scale"
        fi

        # Show output of first iteration.
        my_exec "../globOptVis --show$flag3D --scale $scale --pop-limit $poplimit -p primitives_it$c.bonmin.csv -a $assoc --title \"GlobOpt - $c iteration output\" $visdefparam &"
        
        # Generate relation graphs
        if [ "$doExecPy" = true ]; then
            my_exec "$execPyRelGraph primitives_it$c.bonmin.csv points_primitives_it$c.csv cloud.ply --angles $anglegens --iteration $c"
        fi

        if [ $c -lt $nbExtraIter ]; then
            # Merge adjacent candidates with same dir id. OUT: primitives_merged_it$c.csv, points_primitives_it$c.csv
            my_exec "$executable --merge$flag3D --scale $scale --adopt $adopt --prims primitives_it$c.bonmin.csv -a $assoc --angle-gens $anglegens --patch-pop-limit $poplimit"
        else
            my_exec "../globOptVis --show$flag3D --scale $scale --pop-limit $poplimit -p primitives_it$c.bonmin.csv -a $assoc --title \"GlobOpt - [Dir-Colours] $c iteration output\" $visdefparam --paral-colours --no-rel &"
        fi

        #my_exec "$executable --energy --formulate$flag3D --scale $scale --cloud cloud.ply --unary $unary --pw $pw --cmp 1 --constr-mode $iterationConstrMode --dir-bias $dirbias --patch-pop-limit $poplimit --angle-gens $anglegens --candidates primitives_it$c.bonmin.csv -a $assoc --freq-weight $freqweight  --cost-fn $pwCostFunc $formParams"
    fi
done

if [ "$doExecPy" = true ]; then
    my_exec "$execPyStats .  --angles $anglegens"
fi

#energies
# ../pearl --scale 0.05 --cloud cloud.ply --prims patches.csv --assoc points_primitives.csv --pw 1000 --cmp 1000

# remove safemode
# put back startat
# put back adopt?
