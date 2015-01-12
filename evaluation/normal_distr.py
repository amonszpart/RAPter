
import healpy as hp
import numpy as np

import matplotlib.pyplot as plt

def load_ply(path):
    
    lines = []
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

        
    while (lines[i].split()[0] != 'end_header'):
        i += 1
    
    vstart = i + 1
    
    #read vertices and normals
    for i in range(vstart,vstart+nbV):
        vals = lines[i].split()
        flist = map(float, vals)
        norms.append(flist[0:3]) # load normales from vertex coordinates only
     
    f.close()
    
    return norms

def normalized(a):
    return a / np.linalg.norm(a, 2, 0)



nside = 32

norms = load_ply( "/export/home/kandinsky/nmellado/git/globOpt/data/scenes-paper/bu_lansFull/it30_gaussSphere.ply")

mmap = np.ones((hp.nside2npix(nside)))

for n in norms:
    mmap[hp.vec2pix(nside, n[0], n[1], n[2])] += 1.
#mmap[ hp.vec2pix(nside, norms[0], norms[1], norms[2]) ] += 0.5

mmap /= len(norms)

#print norms[0]

#print hp.vec2pix(nside, norms[0], norms[1], norms[2])

hp.cartview(mmap, nest=True, title="Iteration 30", norm='log', cmap = 'Blues')

plt.savefig('out_it30.svg')
