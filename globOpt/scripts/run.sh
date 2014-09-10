# Executable path
executable="../glob_opt";

###############################################

function print_usage() {
	echo "usage: run.sh scale anglelimit pairwisecost [pop-limit]"
	echo "example: run.sh 0.03 0.2 1"
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

anglegens="90";

visdefparam="" #"--use-tags --no-clusters" #--ids 

echo "angle-limit: $anglelimit"
echo "scale: $scale"
echo "pw: $pw"
echo "pop-limit: $poplimit"

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

# Segmentation. OUTPUT: patches.csv, points_primitives.csv
my_exec "$executable --segment --angle-limit $anglelimit --scale $scale --angle-gens $anglegens"

input="patches.csv";
assoc="points_primitives.csv";


# show segment output
# my_exec "../globOptVis --show --scale $scale --use-tags --ids --pop-limit $poplimit -p patches.csv -a $assoc --title "Segment output"&"
# exit

# Generate candidates. OUT: candidates_it0.csv
my_exec "$executable --generate -sc $scale -al $anglelimit -ald 1 --patch-pop-limit $poplimit -p $input --assoc $assoc --angle-gens $anglegens"
#my_exec "../globOptVis --show --scale $scale --use-tags --ids --pop-limit $poplimit -p candidates_it0.csv -a $assoc --title \"generate output\" &"


# Formulate optimization problem. OUT: "problem" directory
my_exec "$executable --formulate --scale $scale --cloud cloud.ply --unary 10000 --pw $pw --cmp 1 --dir-bias 1 --patch-pop-limit $poplimit --angle-gens $anglegens --candidates candidates_it0.csv -a $assoc"
# Solve optimization problem. OUT: primitives_it0.bonmin.csv
my_exec "$executable --solver bonmin --problem problem -v --time -1 --candidates candidates_it0.csv"

# Show output of first iteration.
my_exec "../globOptVis --show --scale $scale --pop-limit $poplimit -p primitives_it0.bonmin.csv -a $assoc --title \"1st iteration output\" $visdefparam &"

# Merge adjacent candidates with same dir id. OUT: primitives_merged_it0.csv, points_primitives_it0.csv
my_exec "$executable --merge --scale $scale --adopt 0 --prims primitives_it0.bonmin.csv -a $assoc --angle-gens $anglegens"

input="primitives_merged_it0.csv";
assoc="points_primitives_it0.csv";

# Show output of first merge.
my_exec "../globOptVis --show --scale $scale --pop-limit $poplimit -p $input -a $assoc --title \"Merged 1st iteration output\" $visdefparam  &"

# Generate candidates from output of first. OUT: candidates_it1.csv
my_exec "$executable --generate -sc $scale -al 1 -ald 1 --small-mode 2 --patch-pop-limit $poplimit --angle-gens $anglegens -p $input --assoc $assoc"

# Show candidates:
# my_exec "../globOptVis --show --scale $scale -a points_primitives.csv --ids -p candidates_it0.csv --pop-limit $poplimit &"
# Formulate optimization problem. OUT: "problem" directory
my_exec "$executable --formulate --scale $scale --cloud cloud.ply --unary 10000 --pw $pw --cmp 1 --constr-mode 0 --dir-bias 1 --patch-pop-limit $poplimit --angle-gens $anglegens --candidates candidates_it1.csv -a $assoc"

# Solve optimization problem. OUT: primitives_it1.bonmin.csv
my_exec "$executable --solver bonmin -v --problem problem --time -1 --angle-gens $anglegens --candidates candidates_it1.csv"


# Merge adjacent candidates with same dir id. OUT: primitives_merged_it1.csv, points_primitives_it1.csv
my_exec "$executable --merge --scale $scale --adopt 0 --angle-gens $anglegens --prims primitives_it1.bonmin.csv -a $assoc"


# Show output of second iteration.
my_exec "../globOptVis --show --scale $scale --pop-limit 0 --angle-gens $anglegens --prims primitives_merged_it1.csv -a points_primitives_it1.csv --title \"Merged 2nd iteration output\" $visdefparam &"

