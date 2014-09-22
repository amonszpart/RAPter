
# Compute the distance between points of cloud assigned to a primitive
# Return an array with len=len(assign)
def distanceToPrimitives(cloud, assign, primitives):
    return [ [primVar.distanceTo(cloud[a[0]]) for primVar in primitives if primVar.uid == a[1]][0] for a in assign]
