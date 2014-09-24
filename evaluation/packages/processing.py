"""@package Processing
Just a set of various processing function
"""

def removeUnassignedPrimitives(primArray, assignArray):
    
    def isUsed(prim):
        for a in assignArray:
            if a[1] == prim.uid:
                return True
        print "Drop primitive ",prim.uid, prim.did
        return False
        
    return [ x for x in primArray if isUsed(x) ]


def removeUnassignedPoint(primArray, assignArray):
    
    def isUsed(assign):
        for p in primArray:
            if assign[1] == p.uid:
                return True
        print "Drop point ",assign
        return False
        
    return [ x for x in assignArray if isUsed(x) ]
