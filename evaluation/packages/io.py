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
    
    
    
