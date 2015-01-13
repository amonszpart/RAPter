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

 	reprEcho "Running with prims: $rprPrims, assoc: $rprAssoc, iteration: ${rprIter}\n"
	
	rprRepr="representatives_it${rprIter}.csv" # representatives output, contains one patch for each dId
	rprReprAssoc="points_representatives_it${rprIter}.csv" # representatives output, contains associations for representative primitives only
	rprCands="candidates_representatives_it${rprIter}.csv" # candidates generated from representatives
	rprReprOpt="representatives_it${rprIter}.bonmin.csv" # new representatives chosen from candidates
	rprPrimBak="`cutExt $rprPrims`".lvl1.csv


	rprNextId=`expr $c + 1`; #candidates will output here automatically...so we need to know
	rprPw=`my_mult $pw $reprPwMult`
	rprAngLimit=`my_mult $anglelimit 0.5 | awk '{printf "%.3f", $0}'`
	echo "rprAnglimit: $rprAngLimit"

	# representatives
	my_exec "$executable --represent$flag3D -p $rprPrims -a $rprAssoc -sc $scale --cloud cloud.ply --angle-gens $anglegens"
	if ! $dryRun ; then
		echo "mv representatives.csv $rprRepr"
		mv representatives.csv $rprRepr
 		echo "mv points_representatives.csv $rprReprAssoc"
		mv points_representatives.csv $rprReprAssoc
	fi

	
	# ShowRepr
	#cmd="../globOptVis --show$flag3D--scale $scale --pop-limit $poplimit --title \"Representatives\" --angle-gens $angleGens --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-scale --bg-colour .9,.9,.9 --ids --no-rel -p $repr -a $assoc $"
	#my_exec "../globOptVis --show$flag3D --scale $scale --pop-limit $poplimit -p $rprRepr -a $rprReprAssoc --title \"Representatives\" $visdefparam &"

	# Generate from Repr
	if ! $dryRun ; then
		echo "mv candidates_it${rprNextId}.csv candidates_it${rprNextId}_tmp.csv" # move tmp out of the way
		mv candidates_it${rprNextId}.csv candidates_it${rprNextId}_tmp.csv # move tmp out of the way
	fi

	my_exec "$executable --generate$flag3D -sc $scale -al $rprAngLimit -ald ${cand_anglediv} --small-mode 0 --patch-pop-limit $poplimit --angle-gens $candAngleGens --small-thresh-mult $smallThresh -p $rprRepr --assoc $rprReprAssoc --keep-singles"

	if ! $dryRun ; then
		echo "mv candidates_it${rprNextId}.csv $rprCands"
		mv candidates_it${rprNextId}.csv $rprCands
		echo "mv candidates_it${rprNextId}_tmp.csv candidates_it${rprNextId}.csv"
		mv candidates_it${rprNextId}_tmp.csv candidates_it${rprNextId}.csv # move back tmp
	fi
	
	
	# Show candidates
	#my_exec "../globOptVis --show$flag3D --scale $scale --pop-limit $poplimit -p $rprCands -a $rprReprAssoc --title \"GlobOpt-repr_candidates\" $visdefparam &"

	
	# Formulate
	my_exec "$executable --formulate$flag3D $formParams --scale $scale --cloud cloud.ply --unary $unary --pw $rprPw --cmp $cmp --constr-mode patch --dir-bias $dirbias --patch-pop-limit $poplimit --angle-gens $anglegens --candidates $rprCands -a $rprReprAssoc --freq-weight $freqweight --cost-fn $pwCostFunc"

	rprDiagF="diag_it${rprIter}.gv"; 
	rprDiagFTmp="${rprDiagF}RprTmp";
	if ! $dryRun ; then
		echo "cp primitives_it${rprIter}.bonmin.csv primitives_it${rprIter}_rprtmp.csv"
		cp primitives_it${rprIter}.bonmin.csv primitives_it${rprIter}_rprtmp.csv
		if [ -e $rprDiagF ]; then # backup diag_itx.gv
			echo "mv $rprDiagF $rprDiagFTmp";
			mv $rprDiagF "$rprDiagFTmp"
		fi
	fi

	my_exec "$executable --solver$flag3D bonmin --problem problem -v --time -1 --angle-gens $anglegens --bmode $algCode --candidates $rprCands"

	if ! $dryRun ; then
		echo "cp primitives_it${rprIter}.bonmin.csv $rprReprOpt"
		cp primitives_it${rprIter}.bonmin.csv $rprReprOpt
		echo "cp primitives_it${rprIter}_rprtmp.csv primitives_it${rprIter}.bonmin.csv"
		cp primitives_it${rprIter}_rprtmp.csv primitives_it${rprIter}.bonmin.csv
		echo "mv $rprDiagF diag_it${rprIter}.lvl2.gv"
		mv $rprDiagF diag_it${rprIter}.lvl2.gv
		# restore diag_itx.gv
		if [ -e "$rprDiagFTmp" ]; then
			echo "mv $rprDiagFTmp $rprDiagF"
			mv "$rprDiagFTmp" $rprDiagF
		fi
	fi

	#my_exec "../globOptVis --show$flag3D -p $rprReprOpt -a $rprReprAssoc --title \"GlobOpt-RepresentativesOptimized\" --scale $scale --pop-limit $poplimit $visdefparam &"

	# apply representatives - outputs subs.csv
	my_exec "$executable --representBack$flag3D --repr $rprReprOpt -p $rprPrims -a $rprAssoc -sc $scale --cloud cloud.ply --angle-gens $anglegens"

	if ! $dryRun ; then
		echo "mv $rprPrims $rprPrimBak"
		mv $rprPrims $rprPrimBak
		echo "mv subs.csv $rprPrims" #substitue for input
		mv subs.csv $rprPrims
	fi

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