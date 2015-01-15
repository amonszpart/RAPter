
import healpy as hp
import numpy as np
import argparse
import matplotlib.pyplot as plt
from math import sqrt

from deps.plyfile import PlyData, PlyElement

def is_power2(num):
	'states if a number is a power of two'
	return num != 0 and ((num & (num - 1)) == 0)
	
def power2Arg(num):
    if is_power2(num): return num
    msg = "%r is not a power of 2" % num
    raise argparse.ArgumentTypeError(msg)

parser = argparse.ArgumentParser(description='Generate normal distribution from a ply file.')
parser.add_argument('gaussFile')
parser.add_argument('output')
parser.add_argument('title')
parser.add_argument('--nside', default=32, type=power2Arg, help="Power of two used to size the output map.")
parser.add_argument('--useverts', action="store_true", help="Compute the distribution using the vertex coordinates and not the normals (by default).")

args = parser.parse_args()


filename    = args.gaussFile
nside       = args.nside
outfilename = args.output
title       = args.title
useVerts    = args.useverts
print 'Processing ', filename, "(nside=", nside, " useverts=", useVerts , ")"


print "Reading file..."
plydata = PlyData.read(filename)

propX = -1
propY = -1
propZ = -1
elid  = -1

for idx, element in enumerate(plydata.elements):   
    #print element
    if element.name == "vertex":
        elid = idx
        for idp, prop in enumerate(element.properties):
            if(useVerts):
                if prop.name == "x":
                    propX = idp
                if prop.name == "y":
                    propY = idp
                if prop.name == "z":
                    propZ = idp
            else:
                if prop.name == "nx":
                    propX = idp
                if prop.name == "ny":
                    propY = idp
                if prop.name == "nz":
                    propZ = idp
        
if elid == -1 or propX == -1 or propY == -1 or propZ == -1:
    print propX, propY, propZ
    raise RuntimeError("Input field not found")
        
nbVert = len(plydata.elements[elid].data)

print "Building map..."
mmap = np.ones((hp.nside2npix(nside)))
for vertex in plydata.elements[elid].data:
    x = vertex[propX]
    y = vertex[propY]
    z = vertex[propZ]
    norm = sqrt(x*x + y*y + z*z)
    mmap[hp.vec2pix(nside, x/norm, y/norm, z/norm)] += 1.
    

mmap /= nbVert

hp.cartview(mmap, nest=True, title=title, norm='log', cmap = 'Blues')

plt.savefig(outfilename)
