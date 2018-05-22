import numpy as np
from numpy import sqrt, pi, cos, sin, dot, cross, arccos, degrees
from numpy import float as nfl
from numpy.linalg import norm

def __inti__(self):
    pass

def u_vect(vect):
    return(vect/norm(vect))

def angle(A, B, C):
    v1 = B - A
    v2 = B - C
    return(degrees(arccos(np.clip(dot(u_vect(v1), u_vect(v2)), -1.0, 1.0))))

def dih(A, B, C, D):
    r1 = A - B
    r2 = B - C
    r3 = C - D
    n1 = cross(r1, r2)
    n2 = cross(r2, r3)
    return(degrees(arccos(np.clip(dot(u_vect(n1), u_vect(n2)), -1.0, 1.0))))

def geometrical_center(coord):
    return(sum(coord)/float(len(coord)))

def norm_sphere():
    v_sphere = np.random.normal(0.0, 1, (5000,3))
    return(np.array([ u_vect(vect) for vect in v_sphere]))


#####################################
# Bonded Potentials
#####################################

def select_positions(traj,  centers):
    indices = [ atom_indices[i]  for i in centers  ]
    return([np.array(traj.positions[j]) for j in indices])

def SimpleHarmonic(term, traj, centers):
    positions = select_positions(traj, centers)
    dist = np.linalg.norm(positions[0] - positions[1])
    energy = 0.5 * float(term.parameters[1]) * (dist-float(term.parameters[0]))**2.0
    return(energy)

def SimpleHarmonicSquared(term, traj, centers):
    positions = select_positions(traj, atom_indices, centers)
    dist = positions[0] - positions[1]
    energy = 0.5 * float(term.parameters[1]) * (dist**2.0-float(term.parameters[0])**2.0)**2.0
    return(energy)   

def HarmonicAngle(term, traj,  centers):
    positions = select_positions(traj, centers)
    ang = angle(positions[0],positions[1],positions[2])
    print(ang )
    energy = 0.5 * float(term.parameters[1]) * (np.radians(ang)-np.radians(float(term.parameters[0])))**2.0
    return(energy)

def HarmonSqCosAng(term, positions):
    pass

def HarmonDih(term, positions):
    pass

def ReB(term, positions):
    pass

def CBT(term, positions):
    pass

#####################################
# Nonbonded Potentials
#####################################

def LJ(coefs, dist, form):
    if int(form) == 1:
       return(float(coefs[1])/dist**12.0 - float(coefs[0])/dist**6.0)        
    elif int(form) == 2:
       return(4*float(coefs[0])*((float(coef[1])/dist)**12.0-(float(coef[1])/dist)**6.0))
    else:
        raise Exception

def Coul(coefs, dist, eps):
    return(138.935458/(eps)*(float(coefs[0])*float(coefs[1]))/dist**2.0)
