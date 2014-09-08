echo "usage: run.sh scale anglelimit pairwisecost"
scale=$1
echo "scale: $scale"
anglelimit=$2
echo "angle-limit: $anglelimit"
pw=$3
echo "pw: $pw"

poplimit=5

function my_exec() {
	echo $1;
	eval $1
}

executable="../glob_opt";


my_exec "$executable --segment --angle-limit $anglelimit --scale $scale"
my_exec "$executable --generate -sc $scale -al $anglelimit -ald 1 -p patches.csv --patch-pop-limit $poplimit" #-o candidates_it0.csv
my_exec "$executable --formulate --scale $scale --cloud cloud.ply --unary 10000 --pw $pw --cmp 1 --candidates candidates_it0.csv --dir-bias 1 --patch-pop-limit $poplimit"
my_exec "$executable --solver bonmin -v --problem problem --time -1 --candidates candidates_it0.csv"
my_exec "../globOptVis --show --scale $scale -a points_primitives.csv --ids -p primitives_it0.bonmin.csv --pop-limit $poplimit &"

my_exec "$executable --generate -sc $scale -al 1 -ald 1 -p primitives_it0.bonmin.csv --small-mode 2 --patch-pop-limit $poplimit"
my_exec "$executable --formulate --scale $scale --cloud cloud.ply --unary 10000 --pw $pw --cmp 1 --candidates candidates_it1.csv --constr-mode 0 --dir-bias 1  --patch-pop-limit 6"
my_exec "$executable --solver bonmin -v --problem problem --time -1 --candidates candidates_it1.csv"
my_exec "../globOptVis --show --scale $scale -a points_primitives.csv --ids -p primitives_it1.bonmin.csv --pop-limit 0 &"