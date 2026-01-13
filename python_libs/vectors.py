import numpy as np

# distance vector:
def dist(r1,r2):
    r = []
    for r1i,r2i in zip(r1,r2):
        r.append(r2i-r1i)
    return np.array(r)

# absolute distance:
def abs_dist(r1,r2):
    return np.sqrt( np.sum([(r1i-r2i)**2 for r1i,r2i in zip(r1,r2)]) )

# euclidean vector norm:
def norm(r):
    return np.sqrt(np.dot(r,r))

# euclidean vector norm:
def normalize(r):
    return 1/norm(r) * r

# scalar product:
def scalar(u,v):
    return np.dot(u,v)

# cross product
def cross(u,v):
    if len(u)!=3 or len(v)!=3:
        raise RuntimeError("vector length invalid for cross product")
    return np.array([ u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2], u[0]*v[1] - u[1]*v[0] ])

# 3D spherical coordinates vector
def R_sphere(r,phi,costheta):
    sintheta = np.sqrt(1-costheta**2) # sin(x) = sqrt(1-cos^2(x))
    return r * np.array([sintheta*np.cos(phi), sintheta*np.sin(phi), costheta])

# check if two vectors are equal
def isEqual(ra,rb):
    return np.array_equal(ra,rb)
    
# check if two vectors are equal up to given precision
def isRoughlyEqual(ra,rb,decimals):
    return np.all( [np.around(ra[dir],decimals)==np.around(rb[dir],decimals) for dir in range(0,len(ra))] )
    