# Executable path
executable="../glob_opt";

###############################################

function print_usage() {
        echo "usage: run.sh scale anglelimit pairwisecost [pop-limit] [3D]"
        echo "example: run.sh 0.03 0.2 1 25 3D"
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

anglegens="90"; # Desired angle generators in degrees. Default: 90.
nbExtraIter=5;  # iteration count. Default: 2.
dirbias="0";	# not-same-dir-id cost offset. Default: 0. Don't use, if freqweight is on.
freqweight="1"; # dataterm = (freqweight / #instances) * datacost. Default: 0. 1 might be too strong...todo
adopt="0";      # Adopt points argument. Default: 0
# In candidate generation, divide angle limit with this to match copies. Default: 1. Set to 10, if too many candidates (variables).
cand_anglediv="1";# for 3D: "2.5";
# multiply scale by this number to get the segmentation (regionGrowing) spatial distance
segmentScaleMultiplier="1";# for 3D: "2.5";

visdefparam="--use-tags --no-clusters --statuses -1,1" #"--use-tags --no-clusters" #--ids
firstConstrMode="patch" # what to add in the first run formulate. Default: 0 (everyPatchNeedsDirection), experimental: 2 (largePatchesNeedDirection).
iterationConstrMode="patch" # what to add in the second iteration formulate. Default: 0 (everyPatchNeedsDirection), experimental: 2 (largePatchesNeedDirection).
premerge=0
startAt=0
smallThresh="10" # smallThresh * scale is the small threshold
smallThreshlimit="1"


echo "scale: $scale"
echo "pw: $pw"
echo "pop-limit: $poplimit"
echo "anglegens: $anglegens"
echo "angle-limit: $anglelimit"
echo "freq-weight: $freqweight"
echo "region-grow-scale-mult: $segmentScaleMultiplier"
echo "smallThresh: " $smallThresh

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
save_args $0 $@ "--freqweight" $freqweight "--angle-limit" $anglelimit "--segment-scale-mult" $segmentScaleMultiplier "--adopt" $adopt "--dirbias" $dirbias "--cand-anglediv" ${cand_anglediv} "--angle-gens" $anglegens

# show command before run
# stop script when the command fails
function my_exec() {
	echo "__________________________________________________________";
	echo -e "\n\n[CALLING] $1";
    eval $1;
  	if [ "$?" -ne "0" ]; then
	    echo "Error detected ($?). ABORT."
	    exit 1
  	fi
}

function energies() {
        echo -e "\nFirst it energy"
        my_exec "$executable --formulate$flag3D --candidates primitives_it0.bonmin.csv -a points_primitives.csv --energy --scale $scale --cloud cloud.ply --unary 10000 --pw $pw --cmp 1 --dir-bias $dirbias --patch-pop-limit $poplimit --angle-gens $anglegens --constr-mode $firstConstrMode"
        echo -e "\nFirst merged energy"
        my_exec "$executable --formulate$flag3D --candidates primitives_merged_it0.csv -a points_primitives_it0.csv --energy --scale $scale --cloud cloud.ply --unary 10000 --pw $pw --cmp 1 --dir-bias $dirbias --patch-pop-limit $poplimit --angle-gens $anglegens --constr-mode $firstConstrMode"
        echo -e "\nSecond it energy"
        my_exec "$executable --formulate$flag3D --candidates primitives_it1.bonmin.csv -a points_primitives_it0.csv --energy --scale $scale --cloud cloud.ply --unary 10000 --pw $pw --cmp 1 --dir-bias $dirbias --patch-pop-limit $poplimit --angle-gens $anglegens --constr-mode $iterationConstrMode"
        echo -e "\nSecond merge energy"
        my_exec "$executable --formulate$flag3D --candidates primitives_merged_it1.csv -a points_primitives_it1.csv --energy --scale $scale --cloud cloud.ply --unary 10000 --pw $pw --cmp 1 --dir-bias $dirbias --patch-pop-limit $poplimit --angle-gens $anglegens --constr-mode $iterationConstrMode"
}

# [0] Segmentation. OUTPUT: patches.csv, points_primitives.csv
if [ $startAt -eq 0 ]; then
    my_exec "$executable --segment$flag3D --angle-limit $anglelimit --scale $scale --dist-limit-mult $segmentScaleMultiplier --angle-gens $anglegens"
fi

input="patches.csv";
assoc="points_primitives.csv";

# show segment output
my_exec "../globOptVis --show$flag3D --scale $scale --use-tags --pop-limit $poplimit -p patches.csv -a $assoc --normals 5 --title \"Segment output\" --no-clusters &"

# merge before start
if [ $premerge -ne 0 ]; then
    my_exec "$executable --merge$flag3D --scale $scale --adopt $adopt --prims $input -a $assoc --angle-gens $anglegens --patch-pop-limit $poplimit"
    cp patches.csv_merged_it-1.csv $input
    cp points_primitives_it-1.csv $assoc
    my_exec "../globOptVis --show$flag3D --scale $scale --use-tags --pop-limit $poplimit -p patches.csv -a $assoc --normals 5 --title \"PreMerge output\" --no-clusters &"
fi

# Generate candidates. OUT: candidates_it0.csv. #small-mode : small patches don't receive any candidates
my_exec "$executable --generate$flag3D -sc $scale -al $anglelimit -ald ${cand_anglediv} --small-mode 0 --patch-pop-limit $poplimit -p $input --assoc $assoc --angle-gens $anglegens --small-thresh-mult $smallThresh"
#my_exec "../globOptVis --show$flag3D --scale $scale --use-tags --ids --pop-limit $poplimit -p candidates_it0.csv -a $assoc --statuses -1,1 --title \"Generate output\" &"

# Formulate optimization problem. OUT: "problem" directory. constr-mode 2: largePatchesNeedDirectionConstraint
my_exec "$executable --formulate$flag3D --scale $scale --cloud cloud.ply --unary 10000 --pw $pw --cmp 1 --constr-mode $firstConstrMode --dir-bias $dirbias --patch-pop-limit $poplimit --angle-gens $anglegens --candidates candidates_it0.csv -a $assoc --freq-weight $freqweight"

# Solve optimization problem. OUT: primitives_it0.bonmin.csv
my_exec "$executable --solver$flag3D bonmin --problem problem -v --time -1 --candidates candidates_it0.csv"

# Show output of first iteration.
my_exec "../globOptVis --show$flag3D --scale $scale --pop-limit $poplimit -p primitives_it0.bonmin.csv -a $assoc --title \"1st iteration output\" $visdefparam --ids# &"
# !! set flag3D to "3D" and recomment these two lines to work with 3D
#energies
#exit

# Merge adjacent candidates with same dir id. OUT: primitives_merged_it0.csv, points_primitives_it0.csv
my_exec "$executable --merge$flag3D --scale $scale --adopt $adopt --prims primitives_it0.bonmin.csv -a $assoc --angle-gens $anglegens --patch-pop-limit $poplimit"
#cp primitives_it0.bonmin.csv primitives_merged_it0.csv
#cp points_primitives.csv points_primitives_it0.csv

# Show output of first merge.
my_exec "../globOptVis --show$flag3D --scale $scale --pop-limit $poplimit -p primitives_merged_it0.csv -a points_primitives_it0.csv --title \"Merged 1st iteration output\" $visdefparam  &"

for c in $(seq 1 $nbExtraIter)
do
    
    smallThresh=$(($smallThresh / 2))
    
    if [ $smallThresh -lt $smallThreshlimit ]; then
        smallThresh=$smallThreshlimit
        adopt="1"
    fi

    echo "smallThreshMult: " $smallThresh
	  echo "__________________________________________________________";
    echo "Start iteration $c"
    
    prevId=`expr $c - 1`;
    nextId=`expr $c + 1`;

    input="primitives_merged_it$prevId.csv";
    assoc="points_primitives_it$prevId.csv";

    # Generate candidates from output of first. OUT: candidates_it$c.csv. #small-mode 2: small patches receive all candidates
    my_exec "$executable --generate$flag3D -sc $scale -al 1 -ald ${cand_anglediv} --small-mode 0 --patch-pop-limit $poplimit --angle-gens $anglegens -p $input --assoc $assoc --small-thresh-mult $smallThresh"

    # Show candidates:
    # my_exec "../globOptVis --show --scale $scale -a $assoc --ids -p candidates_it$c.csv --pop-limit $poplimit &"
    # Formulate optimization problem. OUT: "problem" directory. --constr-mode 2: largePatchesNeedDirectionConstraint
    my_exec "$executable --formulate$flag3D --scale $scale --cloud cloud.ply --unary 10000 --pw $pw --cmp 1 --constr-mode $iterationConstrMode --dir-bias $dirbias --patch-pop-limit $poplimit --angle-gens $anglegens --candidates candidates_it$c.csv -a $assoc --freq-weight $freqweight"

    # Solve optimization problem. OUT: primitives_it$c.bonmin.csv
    my_exec "$executable --solver$flag3D bonmin -v --problem problem --time -1 --angle-gens $anglegens --candidates candidates_it$c.csv"

    # Show output of first iteration.
    my_exec "../globOptVis --show$flag3D --scale $scale --pop-limit $poplimit -p primitives_it$c.bonmin.csv -a $assoc --title \"$nextId nd iteration output\" $visdefparam &"

    if [ $c -lt $nbExtraIter ]; then
        # Merge adjacent candidates with same dir id. OUT: primitives_merged_it$c.csv, points_primitives_it$c.csv
        my_exec "$executable --merge$flag3D --scale $scale --adopt $adopt --angle-gens $anglegens --prims primitives_it$c.bonmin.csv -a $assoc --patch-pop-limit $poplimit"

        # Show output of second iteration.
        my_exec "../globOptVis --show$flag3D --scale $scale --pop-limit 0 --angle-gens $anglegens --prims primitives_merged_it$c.csv -a points_primitives_it$c.csv --title \"Merged $nextId nd iteration output\" $visdefparam &"
    fi
done

#energies
