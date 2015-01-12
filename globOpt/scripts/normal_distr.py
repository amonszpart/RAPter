
import healpy as hp
import numpy as np

import matplotlib.pyplot as plt

def load_ply(path):
    
    lines = []
    verts = []
    norms = []
        
    
    f = open(path, "r")
    for line in f:
        lines.append(line)
        
    if (lines[0] != "ply\n"):
        return 0
    
    i = 1
    #get number of vertices
    while (lines[i].split()[0] != 'element'):
        i += 1
        
    if (lines[i].split()[1] == 'vertex'):
        nbV = int(lines[i].split()[2])
        print str(nbV) + " vertices"
        
        i += 1
        
        #count number of properties: if 3 no normals, if 6 normals
        nbP = 0
#        while (lines[i].split()[0] == 'property'):
#            nbP += 1
#            i += 1
    
    #if ((lines[i].split()[0] == "element") & (lines[i].split()[1] == "face")):
    #    nbF = int(lines[i].split()[2])
    #    print str(nbF) + " faces"
        
    while (lines[i].split()[0] != 'end_header'):
        i += 1
    
    vstart = i + 1
    
    invertedIndex = [[] for x in xrange(nbV)]
    
    #read vertices and normals
    for i in range(vstart,vstart+nbV):
        vals = lines[i].split()
        flist = map(float, vals)
        verts.append(flist[0:3])
        #if (nbP > 3):
        norms.append(flist[0:3])
        #if (nbP > 6):
        #    curvatures.append(flist[6])
     
    f.close()
    
    return verts, np.swapaxes(norms,0,1)

def normalized(a):
    return a / np.linalg.norm(a, 2, 0)



nside = 32

verts, norms = load_ply( "/export/home/kandinsky/nmellado/git/globOpt/data/scenes-paper/bu_lansFull/cloudRGBNormal_patches_gaussSphere.ply")

mmap = np.zeros((hp.nside2npix(nside)))
mmap[ hp.vec2pix(nside, norms[0], norms[1], norms[2]) ] += 1

mmap = mmap

hp.cartview(mmap, nest=True, title="Mollview image NESTED")
plt.savefig('out_patches.png')
