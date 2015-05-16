#!/usr/bin/python

import os; 
import shutil;
import argparse;

#for root, dirs, files in os.walk("./"):
#    for d in dirs:
# 		print(os.path.join(root, d))

parser = argparse.ArgumentParser();
parser.add_argument("folder", help="Folder to start from")
args = parser.parse_args();
if not args.folder:
	os.exit(1);

def countDids(path):
	counts = {};
	for line in open(path):
		floats = map( float, line.split(",") );
		did = int(floats[7]);
		status = int(floats[8]);
		if ( status == 1 ):
			if ( counts.has_key(did) ):
				counts[did] = counts[did] + 1;
			else:
				counts[did] = 1;
	for key in counts:
		print( "did %d: %d" % (key,counts[key]) );
	return len(counts);

#pattern = "noise4_n0.02"
pattern = args.folder;
dest = "png_%s" % (pattern);
if ( not os.path.exists(dest) ):
	os.mkdir( dest );
pattern = "%s_" % (pattern);
scale = 0.035;
root = "./";
for item in os.listdir(root):
	if os.path.isdir(os.path.join(root,item)):
		if item.find(pattern) != -1:
			os.chdir( item );
			#os.system( "../show.py -s %f -p primitives_it10.bonmin.csv -a points_primitives_it9.csv --2d &" % scale );
			print ( os.system("pwd"));
			primName = "primitives_it10.bonmin.csv";
			if ( os.path.exists("primitives_it10.bonmin.refit.csv") ):
				#os.system("../refit -p primitives_it10.bonmin.csv -a points_primitives_it9.csv --scale 0.035 --target-pop 1000");
				primName = "primitives_it10.bonmin.refit.csv";

			cmd = "../globOptVis --show --scale %f --pop-limit 3 --title \"%s\" --angle-gens 90 --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-scale --bg-colour 1.,1.,1. --no-rel -p %s -a points_primitives_it9.csv --screenshot screenshot.png --vis-size 960,600" % (scale,item,primName);
			print cmd;
			os.system( cmd );
			shutil.copyfile( "screenshot.png", "../%s/%s.png" % (dest,item) );

			key = item[ item.find("pw")+2 : ];
			os.system( "echo %s,%d >> ../%s/stats.csv" % (key,countDids(primName),dest) );

			os.chdir("..");


# ls -t | head -n 40 | grep png
# ls -t | head -n 40 | grep png | xargs cp -r -t /home/bontius/workspace/paper-globOpt-github/trunk/paper/SIGG15/figures/noisePwMatrix/
# cp -r png_noise4_n0.0* -t /home/bontius/workspace/paper-globOpt-github/trunk/paper/SIGG15/figures/noisePwMatrix/