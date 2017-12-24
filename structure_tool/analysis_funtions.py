import random
import math
import argparse
import itertools
import numpy as np
import scipy.optimize as opt
from numpy import sqrt, pi, cos, sin, dot, cross, arccos, degrees
from numpy import float as nfl
from numpy.linalg import norm

def radius_of_gyr(traj):
    N = len(traj)
    diff=np.zeros((N**2))
    count=0
    for i in traj:
        for j in traj:
            diff[count]=dot((i - j),(i-j))
            count = count + 1
    Rg= 1/np.float(N)**2 * sum(diff)
    return(np.float(sqrt(Rg)))

def average_end_end_dist(traj):
    return(0)

def PDI():
    return(0)
