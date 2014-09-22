"""@package IO
Generic input/output functions
"""
import numpy as np

def readPointCloudFromPly(path):
    f = open(path, 'r')
    
    points = []
    
    headerSkipped = False
    for line in f:
        if headerSkipped:
            points.append(np.float32(np.array(line.split(' ')[0:3])))
        else:
            headerSkipped = line.find('end_header') != -1
    
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
    
    
    
