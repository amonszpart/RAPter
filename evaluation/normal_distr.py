
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
parser.add_argument('--noscatter', action="store_true", help="Do not display the distribution as a scatter, useful for noisy data")
parser.add_argument('--gaussFile2', default=None, help="Second input file. Distributions will be merged")

args = parser.parse_args()


filename    = args.gaussFile
filename2   = args.gaussFile2
nside       = int(args.nside)
outfilename = args.output
title       = args.title
useVerts    = args.useverts
clamp_pct   = args.clamppct
no_scatter  = args.noscatter
print 'Processing ', filename, "(nside=", nside, " useverts=", useVerts , ")"


print "Reading file..."
plydata  = PlyData.read(filename)


width  = 400
height = 400
#image = np.array()
projector=hp.projector.CartesianProj(xsize=width,ysize=height)


print "Building map..."
#mmap = np.ones((hp.nside2npix(nside)))


def vec2ij(coords):
    return projector.xy2ij(projector.ang2xy(hp.vec2ang(coords)))

fig = plt.figure(figsize=(float(width)/100., float(height)/100.))
#fig.patch.set_visible(False)


def getScatterData(plydata, prevRun = None):

    class mFunctor(object):
        def __init__(self, scattermap = dict()):
            self.scattermap = scattermap
            
        def process(self, vec3):
            i, j = vec2ij(vec3)
            
            if self.scattermap.get((i,j)) != None:
                self.scattermap[(i,j)] += 1
            else:
                self.scattermap[(i,j)] = 1
        
        def returnValue(self):
            return self.scattermap
    
    if (prevRun == None):
        return processData(plydata, mFunctor())
    else:
        return processData(plydata, mFunctor(prevRun))
        
        


def getZImageData(plydata, zImage):

    class mFunctor(object):            
        def process(self, vec3):
            i, j = vec2ij(vec3)
            zImage[i,j] += 1.
        
        def returnValue(self):
            return zImage
    

    return processData(plydata, mFunctor())

def processData(plydata, functor):
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
        
    for vertex in plydata.elements[elid].data:
        x = vertex[propX]
        y = vertex[propY]
        z = vertex[propZ]
        norm = sqrt(x*x + y*y + z*z)        
        functor.process (np.array([x/norm, y/norm, z/norm]))

    return functor.returnValue()
    
print "Processing point cloud..."    
    
if no_scatter:
    zImage = np.empty([height, width])
    zImage.fill(0.0001)
    
    getZImageData(plydata, zImage) 
    
    if filename2 != None: 
        print "Processing second point cloud..."  
        getZImageData(PlyData.read(filename2), zImage)
       
    plt.imshow(zImage, cmap='Blues', norm=LogNorm(vmin=np.min(zImage), vmax=np.max(zImage)))
else:
    scattermap = getScatterData (plydata)
    
    
    if filename2 != None: 
        print "Processing second point cloud..."  
        scattermap = getScatterData(PlyData.read(filename2), scattermap)
    
    print "Building scatter..."   
    ## convert scattermap to scatter arrays
    sx = []
    sy = []
    sz = []
    for key, value in scattermap.iteritems():
        sx.append(key[1])
        sy.append(height-key[0]-1)
        sz.append(value)
    #plt.scatter(sx, sy, s=sz, alpha=0.5)
    plt.scatter(sx, sy, s=40, c=sz, alpha=0.5, cmap='Blues')

ax1 = plt.gca()
ax1.set_xlim(0., width)
ax1.set_ylim(0., height)
ax1.axis('off')

plt.savefig(outfilename, dpi=300, bbox_inches='tight' )

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
