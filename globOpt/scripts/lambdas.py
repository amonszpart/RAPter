#!/usr/bin/python3
import numpy as np
import math

def fFrom(val):
	return 3+math.log10(val);

a = np.float64(1e5)
#bs = [np.float64(0.1), np.float64(10), np.float64(100), np.float64(1e5), np.float64(1e8)]
pws = [ "1e-3", "1e-2", "1e-1", "1e0", "1e1", "2.5e1", "5e1", "7.5e1", "1e2", "2.5e2", "5e2", "7.5e2", "1e3", "5e3", "1e4", "5e4", "1e5", "5e5", "1e6", "5e6", "1e7", "5e7", "1e8", "1e10" ];
pws = [ "5e-3", "0.1","10","100","5e3", "1e5","5e6","1e8"];
bs = map( float, pws );
ab = []
div = np.float64(1e6)

for b in bs:
    #ab.append( [(a/(a+b), b/(a+b))] );
    a2 = fFrom(a)
    b2 = fFrom(b)
    ab.append([a2 / (a2 + b2), b2 / (a2 + b2)])

print("pw0\t", end='')
for v in bs:
    print(("%.1f & " % (v)).strip("0"), end='')
print("")

#print("data:\t", end='')
for pair in ab:
    print("%.2f & " % (pair[0]), end='')
print("")


#print("pw:\t", end='')
for pair in ab:
    print((	"%.2f & " % (pair[1])), end='')
print("")
