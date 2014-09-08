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

echo "angle-limit: $anglelimit"
echo "scale: $scale"
echo "pw: $pw"
echo "pop-limit: $poplimit"

# show command before run
function my_exec() {
	echo $1;
	eval $1
}

# Segmentation. OUTPUT: patches.csv, points_primitives.csv
my_exec "$executable --segment --angle-limit $anglelimit --scale $scale"
# Generate candidates. OUT: candidates_it0.csv
my_exec "$executable --generate -sc $scale -al $anglelimit -ald 1 -p patches.csv --patch-pop-limit $poplimit"
# Formulate optimization problem. OUT: "problem" directory
my_exec "$executable --formulate --scale $scale --cloud cloud.ply --unary 10000 --pw $pw --cmp 1 --candidates candidates_it0.csv --dir-bias 1 --patch-pop-limit $poplimit"
# Solve optimization problem. OUT: primitives_it0.bonmin.csv
my_exec "$executable --solver bonmin -v --problem problem --time -1 --candidates candidates_it0.csv"
# Show output of first iteration.
my_exec "../globOptVis --show --scale $scale -a points_primitives.csv --ids -p primitives_it0.bonmin.csv --pop-limit $poplimit &"

# Generate candidates from output of first. OUT: candidates_it1.csv
my_exec "$executable --generate -sc $scale -al 1 -ald 1 -p primitives_it0.bonmin.csv --small-mode 2 --patch-pop-limit $poplimit"
# Show candidates:
# my_exec "../globOptVis --show --scale $scale -a points_primitives.csv --ids -p candidates_it0.csv --pop-limit $poplimit &"
# Formulate optimization problem. OUT: "problem" directory
my_exec "$executable --formulate --scale $scale --cloud cloud.ply --unary 10000 --pw $pw --cmp 1 --candidates candidates_it1.csv --constr-mode 0 --dir-bias 1  --patch-pop-limit 6"
# Solve optimization problem. OUT: primitives_it1.bonmin.csv
my_exec "$executable --solver bonmin -v --problem problem --time -1 --candidates candidates_it1.csv"
# Show output of second iteration.
my_exec "../globOptVis --show --scale $scale -a points_primitives.csv --ids -p primitives_it1.bonmin.csv --pop-limit 0 &"