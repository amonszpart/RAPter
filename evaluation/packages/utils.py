
# Compute the distance between points of cloud assigned to a primitive
# Return an array with len=len(assign)
def distanceToPrimitives(cloud, assign, primitives):
    return [ [primVar.distanceTo(cloud[a[0]]) for primVar in primitives if primVar.uid == a[1]] for a in assign]
    
    
   
import packages.orderedSet as orderedSet 
def parseAngles(strAngle):
    angles = orderedSet.OrderedSet()
    angles.add(0.)
    if len(strAngle) == 1:
        strAngle = strAngle[0].split(',')

    for genAngle in strAngle:
        a = float(genAngle)
        while a <= 180.:
            angles.add(a)
            a+= float(genAngle)
    
    return angles
    
