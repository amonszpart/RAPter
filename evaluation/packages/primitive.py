"""@package Primitive
This module is the Python counterpart of the C++ LinePrimitive and and PlanePrimitive classes

See C++ GlobOpt project for more details on Primitives
"""
import numpy as np

class Primitive(object):

    def __init__(self, 
                 uid,
                 did,
                 pos=np.zeros(3), 
                 normal=np.array([ 1.,  0.,  0.]) ):
        self.uid    = uid
        self.did    = did
        self.pos    = pos
        self.normal = normal
    
    def distanceTo(self, pos) :
        return np.dot(pos-self.pos, self.normal)
        
def readPrimitivesFromFile(path):
    f = open(path, 'r')
    
    primitives = []
    
    for line in f:
        if line[0] != '#':
            seqline = line.split(',')
            
            pos    = np.float32(np.array(seqline[0:3]))
            normal = np.float32(np.array(seqline[3:6]))
            uid    = int(seqline[6])
            did    = int(seqline[7])
            
            primitives.append(Primitive( uid, did, pos, normal ))

    return primitives
