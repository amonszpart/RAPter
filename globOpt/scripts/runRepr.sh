#! /bin/bash

function reprEcho()
{
	echo -e "[runRepr] $1\n" >&1;
}

function printReprUsage()
{
	reprEcho "\n\nrunRepr usage: primitives_it\?.bonmin.csv points_primitives_it\?-1.csv iteration\n";
	exit 1
}

function runRepr() {
	if [[ -z "$1" ]]; then printReprUsage; exit 1; else rprPrims=$1; fi
	if [[ -z "$2" ]]; then printReprUsage; exit 1; else rprAssoc=$2; fi
	if [[ -z "$3" ]]; then printReprUsage; exit 1; else rprIter=$3; fi	

 	reprEcho "Running with prims: $rprPrims, assoc: $rprAssoc, iteration: $rprIter\n"
	
	rprRepr="representatives_it$rprIter.csv" # representatives output, contains one patch for each dId
	rprReprAssoc="points_representatives_it$rprIter.csv" # representatives output, contains associations for representative primitives only
	rprCands="candidates_representatives_it$rprIter.csv" # candidates generated from representatives
	rprReprOpt="representatives_it$rprIter.bonmin.csv" # new representatives chosen from candidates
	rprPrimBak="`cutExt $rprPrims`".lvl1.csv

	#prevId=`expr $c - 1`;
	rprNextId=`expr $c + 1`; #candidates will output here automatically...so we need to know
	rprPw=`my_mult $pw $reprPwMult`

	# representatives
	#cmd="../glob_opt --represent$flag3D -p $prims -a $assoc -sc $scale --cloud cloud.ply --angle-gens $angleGens"
	my_exec "$executable --represent$flag3D -p $rprPrims -a $rprAssoc -sc $scale --cloud cloud.ply --angle-gens $anglegens"
	mv representatives.csv $rprRepr
	mv points_representatives.csv $rprReprAssoc

	# ShowRepr
	#cmd="../globOptVis --show$flag3D--scale $scale --pop-limit $poplimit --title \"Representatives\" --angle-gens $angleGens --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-scale --bg-colour .9,.9,.9 --ids --no-rel -p $repr -a $assoc $"
	#my_exec "../globOptVis --show$flag3D --scale $scale --pop-limit $poplimit -p $rprRepr -a $rprReprAssoc --title \"Representatives\" $visdefparam &"

	# Generate from Repr
	#cmd="../glob_opt --generate$flag3D -sc $scale -al $anglelimit -ald 1 --small-mode 0 --patch-pop-limit $poplimit --angle-gens $angleGens --small-thresh-mult 1 -p $repr --assoc $assoc"
	mv candidates_it$rprNextId.csv candidates_it$rprNextId_tmp.csv # move tmp out of the way
	my_exec "$executable --generate$flag3D -sc $scale -al $anglelimit -ald ${cand_anglediv} --small-mode 0 --patch-pop-limit $poplimit --angle-gens $anglegens --small-thresh-mult $smallThresh -p $rprRepr --assoc $rprReprAssoc --keep-singles"
	echo "mv candidates_it$rprNextId.csv $rprCands"
	mv candidates_it$rprNextId.csv $rprCands
	echo "mv candidates_it$rprNextId_tmp.csv candidates_it$rprNextId.csv"
	mv candidates_it$rprNextId_tmp.csv candidates_it$rprNextId.csv # move back tmp
	
	# Show candidates
	#my_exec "../globOptVis --show$flag3D --scale $scale --pop-limit $poplimit -p $rprCands -a $rprReprAssoc --title \"GlobOpt-repr_candidates\" $visdefparam &"

	# Formulate
	my_exec "$executable --formulate$flag3D $formParams --scale $scale --cloud cloud.ply --unary $unary --pw $rprPw --cmp $cmp --constr-mode patch --dir-bias $dirbias --patch-pop-limit $poplimit --angle-gens $angleGens --candidates $rprCands -a $rprReprAssoc --freq-weight $freqweight --cost-fn $pwCostFunc"

	cp primitives_it$rprIter.bonmin.csv primitives_it$rprIter_rprtmp.csv
	my_exec "$executable --solver$flag3D bonmin --problem problem -v --time -1 --angle-gens $anglegens --bmode $algCode --candidates $rprCands"
	cp primitives_it$rprIter.bonmin.csv $rprReprOpt
	cp primitives_it$rprIter_rprtmp.csv primitives_it$rprIter.bonmin.csv

	#my_exec "../globOptVis --show$flag3D -p $rprReprOpt -a $rprReprAssoc --title \"GlobOpt-RepresentativesOptimized\" --scale $scale --pop-limit $poplimit $visdefparam &"

	# apply representatives - outputs subs.csv
	my_exec "$executable --representBack$flag3D --repr $rprReprOpt -p $rprPrims -a $rprAssoc -sc $scale --cloud cloud.ply --angle-gens $anglegens"
	mv $rprPrims $rprPrimBak
	mv subs.csv $rprPrims #substitue for input

	#my_exec "../globOptVis --show$flag3D --scale $scale --pop-limit $poplimit --title \"GlobOpt - ReprBack output\" --angle-gens $angleGens --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-scale --bg-colour .3,.3,.3 --ids --no-rel -p subs.csv -a $oldAssoc"
	#my_exec "../globOptVis --show$flag3D -p $rprPrims -a $rprAssoc --title \"GlobOpt-ReprOutput\" --scale $scale --pop-limit $poplimit $visdefparam &"
	
	# ../glob_opt --energy --formulate --scale 0.015 --cloud cloud.ply --unary 10000 --pw 15 --cmp $cmp --constr-mode patch --dir-bias 0 --patch-pop-limit 5 --angle-gens 60,90 --candidates primitives_it.bonmin.csv -a --freq-weight 1000 --cost-fn spatsqrt

	# showmy
	# ../globOptVis --show --scale 0.02 --pop-limit 5 --title "GlobOpt - ReprOpt output" --angle-gens 60,90 --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-scale --bg-colour .3,.3,.3 --ids --no-rel -p representatives.bonmin.csv.my -a points_representatives.csv --show-spatial
}





	#out=representatives.bonmin.csv
	#cmd="../glob_opt --energy --formulate$flag3D --scale $scale --cloud cloud.ply --unary $unary --pw $pw --cmp $cmp --constr-mode patch --dir-bias 0 --patch-pop-limit $poplimit --angle-gens $angleGens --candidates $out -a $assoc --freq-weight $freqweight --cost-fn spatsqrt"

	#out=representatives.bonmin.csv.my
	#cmd="../glob_opt --energy --formulate$flag3D --scale $scale --cloud cloud.ply --unary $unary --pw $pw --cmp $cmp --constr-mode patch --dir-bias 0 --patch-pop-limit $poplimit --angle-gens $angleGens --candidates $out -a $assoc --freq-weight $freqweight --cost-fn spatsqrt"
	#echo $cmd
	#`$cmd >log.tmp`
	#cat log.tmp