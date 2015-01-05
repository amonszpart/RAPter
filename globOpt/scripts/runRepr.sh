#! /bin/bash

## These are scripts how to run a representative problem

scale=0.02
unary=1000
pw=10000
poplimit=5
angleGens="60,90"
cands="candidates_it0.csv"
pwCostFunc="spatsqrt"
freqweight=1000
flag3D=""
prims="primitives_it6.bonmin.csv"
assoc="points_primitives_it5.csv"
anglelimit=0.1
algCode=0

# representatives
cmd="../glob_opt --represent -p $prims -a $assoc -sc $scale --cloud cloud.ply"
echo $cmd
`$cmd >/dev/null 2>/dev/null`
repr="representatives.csv"

# ShowRepr
cmd="../globOptVis --show --scale $scale --pop-limit $poplimit --title \"Representatives\" --angle-gens $angleGens --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-scale --bg-colour .9,.9,.9 --ids --no-rel -p $repr -a $assoc"
echo $cmd
`$cmd >/dev/null 2>/dev/null`

# Generate from Repr
cmd="../glob_opt --generate -sc $scale -al $anglelimit -ald 1 --small-mode 0 --patch-pop-limit $poplimit --angle-gens $angleGens --small-thresh-mult 1 -p $repr --assoc $assoc"
echo $cmd
`$cmd >/dev/null 2>/dev/null`

# Show candidates
cp candidates_it0.csv candidates_it0_tmp.csv
cmd="../globOptVis --show --scale $scale --pop-limit $poplimit -p candidates_it0.csv -a points_representatives.csv --title \"GlobOpt-repr_candidates\" --angle-gens $angleGens --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-scale --bg-colour .9,.9,.9 --ids --no-rel"
echo $cmd
`$cmd >/dev/null 2>/dev/null`
cp candidates_it0.csv candidates_representatives.csv
cp candidates_it0_tmp.csv candidates_it0.csv
cands="candidates_representatives.csv"
assoc="points_representatives.csv"

# Formulate
cmd="../glob_opt --formulate --scale $scale --cloud cloud.ply --unary $unary --pw $pw --cmp 1 --constr-mode patch --dir-bias 0 --patch-pop-limit $poplimit --angle-gens $angleGens --candidates $cands -a $assoc --freq-weight $freqweight --cost-fn $pwCostFunc"
echo $cmd
`$cmd >/dev/null 2>/dev/null`

cp primitives_it0.bonmin.csv primitives_it0_tmp.csv
cmd="../glob_opt --solver$flag3D bonmin --problem problem -v --time -1 --angle-gens $angleGens --bmode $algCode --candidates $cands"
echo $cmd
`$cmd >log.tmp 2>/dev/null`
cp primitives_it0.bonmin.csv representatives.bonmin.csv
cp primitives_it0_tmp.csv primitives_it0.bonmin.csv
cat log.tmp

cmd="../globOptVis --show --scale $scale --pop-limit $poplimit --title \"GlobOpt - ReprOpt output\" --angle-gens $angleGens --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-scale --bg-colour .3,.3,.3 --ids --no-rel -p representatives.bonmin.csv -a points_representatives.csv --show-spatial"
echo $cmd
`$cmd >/dev/null 2>/dev/null`

out=representatives.bonmin.csv
cmd="../glob_opt --energy --formulate --scale $scale --cloud cloud.ply --unary $unary --pw $pw --cmp 1 --constr-mode patch --dir-bias 0 --patch-pop-limit $poplimit --angle-gens $angleGens --candidates $out -a $assoc --freq-weight $freqweight --cost-fn spatsqrt"
echo $cmd
`$cmd >log.tmp`
cat log.tmp

out=representatives.bonmin.csv.my
cmd="../glob_opt --energy --formulate --scale $scale --cloud cloud.ply --unary $unary --pw $pw --cmp 1 --constr-mode patch --dir-bias 0 --patch-pop-limit $poplimit --angle-gens $angleGens --candidates $out -a $assoc --freq-weight $freqweight --cost-fn spatsqrt"
echo $cmd
`$cmd >log.tmp`
cat log.tmp

# ../glob_opt --energy --formulate --scale 0.015 --cloud cloud.ply --unary 10000 --pw 15 --cmp 1 --constr-mode patch --dir-bias 0 --patch-pop-limit 5 --angle-gens 60,90 --candidates primitives_it.bonmin.csv -a --freq-weight 1000 --cost-fn spatsqrt

# showmy
# ../globOptVis --show --scale 0.02 --pop-limit 5 --title "GlobOpt - ReprOpt output" --angle-gens 60,90 --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-scale --bg-colour .3,.3,.3 --ids --no-rel -p representatives.bonmin.csv.my -a points_representatives.csv --show-spatial
