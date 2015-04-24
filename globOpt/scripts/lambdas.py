#!/usr/bin/python

a = 1e5;
bs = [ 0.1, 10, 100, 1e5, 1e8 ];

for b in bs:
	print( "%f %f" % (a/(a+b), b/(a+b)) );