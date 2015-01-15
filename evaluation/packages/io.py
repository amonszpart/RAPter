"""@package IO
Generic input/output functions
"""
import numpy as np
from deps.plyfile import PlyData, PlyElement

def readPointCloudFromPly(path):
    plydata = PlyData.read(path)
    
    propX = -1
    propY = -1
    propZ = -1
    elid  = -1

    for idx, element in enumerate(plydata.elements):
        if element.name == "vertex":
            elid = idx
            for idp, prop in enumerate(element.properties):
                if prop.name == "x":
                    propX = idp
                if prop.name == "y":
                    propY = idp
                if prop.name == "z":
                    propZ = idp
                    
    if elid == -1 or propX == -1 or propY == -1 or propZ == -1:
        raise RuntimeError("IO: Input field not found")
    
    points = []
    for vertex in plydata.elements[elid].data:
        points.append(np.float32(np.array([vertex[propX], vertex[propY], vertex[propZ]])))
    
    return points
    
def readPointAssignementFromFiles(path):
    f = open(path, 'r')
    
    assignement = []
    
    for line in f:
        if line[0] != '#':
            assignement.append(np.int16(line.split(',')[0:2]))
    
    return assignement
    
def readPrimitiveCorrespondancesFromFiles(path, primset1, primset2):
    f = open(path, 'r')
    
    corresp    = {}
    correspUid = {}
    
    for line in f:
        if line[0] != '#':
            fline = np.int16(line.split(',')[0:6])
            
            p1 = None
            p2 = None
            
            for p in primset1:
                if p.uid == fline[0] and p.did == fline[2]:
                    p1 = p
                    break
            
            for p in primset2:
                if p.uid == fline[3] and p.did == fline[5]:
                    p2 = p
                    break
            
            if (p1 != None and p2 != None):
                corresp[p1] = p2
                correspUid[p1.uid] = p2.uid
            else: 
                print "Cannot find ",fline
            
    f
    
    return corresp, correspUid
    
    
    
