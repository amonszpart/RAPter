
import healpy as hp
import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from math import sqrt

from deps.plyfile import PlyData, PlyElement

def is_power2(num):
    'states if a number is a power of two'
    return num == 2 or ( num%2==0 and is_power2(num/2))
	

	
def power2Arg(num):
    if int(num) != 0:
        if is_power2(int(num)): return num
    msg = "%r is not a power of 2" % num
    raise argparse.ArgumentTypeError(msg)

parser = argparse.ArgumentParser(description='Generate normal distribution from a ply file.')
parser.add_argument('gaussFile')
parser.add_argument('output')
parser.add_argument('title')
parser.add_argument('--nside', default=32, type=power2Arg, help="Power of two used to size the output map.")
parser.add_argument('--clamppct', default=-1., type=float, help="Max percentage used to clamp the distribution.")
parser.add_argument('--useverts', action="store_true", help="Compute the distribution using the vertex coordinates and not the normals (by default).")

args = parser.parse_args()


filename    = args.gaussFile
nside       = int(args.nside)
outfilename = args.output
title       = args.title
useVerts    = args.useverts
clamp_pct   = (args.clamppct)
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

print hp.nside2npix(nside), hp.nside2npix(nside) / 2, sqrt(hp.nside2npix(nside) / 2)
h=sqrt(hp.nside2npix(nside) / 2)


#image = np.array()
projector=hp.projector.CartesianProj(xsize=800,ysize=400)


print "Building map..."
mmap = np.ones((hp.nside2npix(nside)))

scattermap = dict()

for vertex in plydata.elements[elid].data:
    x = vertex[propX]
    y = vertex[propY]
    z = vertex[propZ]
    norm = sqrt(x*x + y*y + z*z)
    
    i, j =  projector.xy2ij(projector.ang2xy(hp.vec2ang(np.array([x/norm, y/norm, z/norm]))))
    
    if scattermap.get((i,j)) != None:
        scattermap[(i,j)] += 1
    else:
        scattermap[(i,j)] = 1
    
    

    #print i, j
    #exit()
    #mmap[hp.vec2pix(nside, x/norm, y/norm, z/norm)] += 1.
print "Building scatter..."
    
## convert scattermap to scatter arrays
sx = []
sy = []
sz = []
for key, value in scattermap.iteritems():
    sx.append(key[0])
    sy.append(key[1])
    sz.append(value)
    
  
plt.figure(figsize=(5.,2.5))
#plt.scatter(sx, sy, s=sz, alpha=0.5)
plt.scatter(sx, sy, s=40, c=sz, alpha=0.5, cmap='Blues')

ax1 = plt.gca()
ax1.get_xaxis().set_ticks([])
ax1.get_yaxis().set_ticks([])

plt.savefig(outfilename, dpi=300, bbox_inches='tight', )

exit()


print len(scattermap)

mmap /= float(nbVert) / 100.

mini = np.min(mmap)
maxi = np.max(mmap)

print "Distribution max: ", maxi

if clamp_pct > 0.: 
    clamp_pct *= 100.
    
    if clamp_pct >= maxi:
        print "Clamp is bigger than distribution max. Ignore"
    else:
        mmap = np.clip(mmap, mini, clamp_pct)
        print "Clamp to: ", clamp_pct
    maxi = clamp_pct   
    

hp.cartview(mmap, nest=True, title=title, cmap = 'Blues', unit='%')

plt.savefig(outfilename, dpi=300)
