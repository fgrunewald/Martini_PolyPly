import random
import math
import argparse
import itertools
import numpy as np
import scipy.optimize as opt
from numpy import sqrt, pi, cos, sin, dot, cross, arccos, degrees
from numpy import float as nfl
from numpy.linalg import norm
from Martini_PolyPly.structure_tool.mc_poly_growth import *
from Martini_PolyPly.structure_tool.analysis_funtions import *
from Martini_PolyPly.structure_tool.geometrical_functions import *
from Martini_PolyPly.structure_tool.force_field_tools import *

def accaptable(E, temp, prev_E, overlap_limit):
    if E < prev_E:
       return(True)
    elif E > overlap_limit:
       return(False)
    else:
       N = np.random.normal(0,1,(1,1))
       F = np.exp(- (E-prev_E)/(kb*temp))
       if N < F:
          return(True)
       else:
          return(False)

def determine_step_legnth(coords, bb_indices):
       bb_coord = [coords[i] for i in bb_indices]
       step_length =  norm(sum(bb_coord))
       g0 = geometrical_center(coords)
       max_dist_from_center = max([ norm(g0 - point) for point in coords])
       size =  max_dist_from_center + 0.43
       nexcl = math.ceil(size/step_length)
       print("step:",step_length)
       print("size:", size)
       print("Will exclude",nexcl, "interactions.")
       return(step_length, size, nexcl)

def take_step(vectors, step_length, item):
    index = random.randint(0, len(vectors) - 1)
    new_item = item + vectors[index] * step_length
    return(new_item, index)

def Hamiltonion(ff, traj, display):
    traj = traj.reshape(-1,3)
    bonded = bonded_pot(ff, traj)
    angle = angle_pot(ff, traj)
    dihedral = dihedral_pot(ff, traj)
    vdw = Vdw_pot(ff, traj)
    if display:
       for term, name in zip([bonded, angle, dihedral, vdw],['bonds', 'angle', 'dihedral', 'vdw']):
           print(name, term)
    return(bonded + angle + dihedral + vdw)


def metropolis_monte_carlo(ff, conf, temp, n_repeat, step_length, max_steps, verbose):
    traj = np.array(conf)
    count = 0
    prev_E = 0.0
    rejected = 0
    print('\n++++++++++ STARTING MONTE CARLO MODULE +++++++++\n')
    while count < n_repeat -1:
       #print('----------------')     
       while True:
          vector_bundel = norm_sphere()
          new_coord, index = take_step(vector_bundel, step_length, traj[len(traj) - 1])
          new_traj = np.append(traj, np.array([new_coord])).reshape(-1,3)
          total_E  = Hamiltonion(ff, new_traj, verbose)/len(new_traj)
          atom_type =  ff['atoms'][count]['typ']
          epsilon = ff['nonbond_params'][(atom_type, atom_type)]['epsilon']
          sigma = ff['nonbond_params'][(atom_type, atom_type)]['sigma']
          overlap_limit = LJ(sigma , epsilon, sigma * 0.8)      
          if accaptable(total_E, temp, prev_E, overlap_limit):
            if verbose:
               print('accapted')
               print(total_E * len(new_traj))
            prev_E = total_E
            traj = new_traj
            print(count)
            count = count + 1
            break
          elif count < max_steps:
            #print('rejected')
            if verbose:
               print('rejected')         
            rejected = rejected + 1
            vector_bundel = np.delete(vector_bundel, index, axis=0)
          else:
            print('+++++++++++++ FATAL ERROR ++++++++++++++++++++++++++++\n')
            print('Exceeded maximum number of steps in the monte-carlo module.')
            print('If you know what you do set -maxconstraints to -1')

    print('++++++++++++++++ RESULTS FORM MONTE CARLO MODULE ++++++++++++++\n')
    print('Total Energy:', Hamiltonion(ff, new_traj, verbose))
    print('Radius of Gyration:', radius_of_gyr(traj))
    print('Number of rejected MC steps:', rejected)       
    return(traj)

def generate_chains(ff, conf, step_length, nexcl, chain, n_mon, temp, max_steps, conv, verbose):
    traj = conf 
    if chain in '[ linear-random ]':
       traj = metropolis_monte_carlo(ff, conf, temp, n_mon, step_length, max_steps, verbose)
    elif chain in '[ PEGylated-bilayer ]':
       traj = metropolis_monte_carlo(ff,  conf,temp, n_mon, step_length, max_steps, verbose)
    elif chain in '[ Polystyrene ]':
       print('PS module')
    return(traj)

def add_environment():
    return(None)

def build_system(topfile, structure_file, chain, env, n_chains, n_mon, box_vect, temp, max_steps, conv, verbose):
    ff = read_itp(topfile)
    # Since FF is global there is no output of convert_constraints  
    verbose=False
    ff = convert_constraints(ff, conv)
    conf = read_conf_file(structure_file, 'gro')
    step_length, size, nexcl = determine_step_legnth(conf, [0])
    traj = generate_chains(ff, conf, step_length, nexcl, chain, n_mon, temp, max_steps, conv, verbose)
    #print(traj)
    write_gro_file(traj,'out.gro',len(traj))
    return(None)   
