#!/usr/bin/python
import argparse;
import numpy as np;
import math;

parser = argparse.ArgumentParser();
parser.add_argument("file", help="File to refit")
args = parser.parse_args();

in2cluster = {};
clusters = [ ];
limitRad = 0.01 / 180.0 * math.pi;
lineId = 0;
lines = [];
for line in open(args.file):
	lines.append(line);
	floats = map( float, line.split(",") );
	n = np.array( [ floats[3], floats[4], floats[5] ] );
	print "\nnew:",n;
	hit = False;
	for i,cluster in enumerate(clusters):
		for normal in cluster:
			print "normal: ", normal
			try:
				ddot = n.dot( normal );
				angRad = math.acos( ddot );
				#print n.cross(normal).norm()
				#np.cross(n,normal)
				#angRad = math.atan2( np.cross(n,normal).norm(), n.dot(normal) );
				break
			except:
				angRad = 0.;

			if ( abs(angRad) < limitRad ):
				print "ok", n, normal
				hit = True;
				in2cluster[ lineId ] = i;
				if ddot < 0.:
					n = -1. * n;
				break;
			else:
				print "no: ", angRad, " > ", limitRad, "diff: ", angRad - limitRad
	if not hit:
		print "enqueing: ", n;
		clusters.append( [n] );
		in2cluster[ lineId ] = len(clusters)-1;
	else:
		clusters[ in2cluster[lineId] ].append( n );

	lineId += 1;

avgs = [];
for cluster in clusters:
	avg = np.array( [0.,0.,0.] );
	for n in cluster:
		avg += n;
	avgs.append( avg / float(len(cluster)) );

f2 = open("out.csv",'w');
for key in in2cluster:
	outNormal = avgs[ in2cluster[key] ];
	print key, "->", in2cluster[key], ": ", lines[key], "->", outNormal;
	arr = lines[key].split(",");
	line2 = "%s,%s,%s,%f,%f,%f" % ( arr[0], arr[1], arr[2], outNormal[0], outNormal[1], outNormal[2] );
	for i in range(6,len(arr)):
		line2 = "%s,%s" % (line2,arr[i]);
	f2.write( line2 );
f2.close()
