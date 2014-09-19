"""@package Processing
Just a set of various processing function
"""

def removeUnassignedPrimitives(primArray, assignArray):
    
    def isUsed(prim):
        for a in assignArray:
            if a[1] == prim.uid:
                return True
        return False
        
    return [ x for x in primArray if isUsed(x) ]


def removeUnassignedPoint(primArray, assignArray):
    
    def isUsed(assign):
        for p in primArray:
            if assign[1] == p.uid:
                return True
        return False
        
    return [ x for x in assignArray if isUsed(x) ]
