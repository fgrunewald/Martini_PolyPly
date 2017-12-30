import random
import math
import argparse
import itertools
import numpy as np
import scipy.optimize as opt
import multiprocessing
from numpy import sqrt, pi, cos, sin, dot, cross, arccos, degrees
from numpy import float as nfl
from numpy.linalg import norm
from Martini_PolyPly.structure_tool.mc_poly_growth import *
from Martini_PolyPly.structure_tool.analysis_funtions import *
from Martini_PolyPly.structure_tool.geometrical_functions import *
from Martini_PolyPly.structure_tool.force_field_tools import *
from Martini_PolyPly.structure_tool.environment import *


def accaptable(E, temp, prev_E):
    if E == math.inf:
       return(False)
    if E < prev_E:
       return(True)
    else:
       N = np.random.normal(0,1,(1,1))
       F = np.exp(- (E-prev_E)/(kb*temp))
       if N < F:
          return(True)
       else:
          return(False)

def remove_overlap(new_point, traj, tol, sol):
    count=0
    for index, coords in enumerate(traj):
        for point in coords:
       #   print(norm(point-new_point))
          if norm(point - new_point) < tol:
             del traj[index]
             count = count + 1
    print('Removed', count, 'atoms')
    return(traj)

def is_overlap(new_point, traj, tol, nexcl=1):
    n = len(traj) - nexcl
    distances = [ norm(point - new_point) for point in traj[0:n]]
    return( any([ dist <  tol for dist in distances]))


def constraints(new_point, list_of_constraints):
    status = []
   
    for const in list_of_constraints:
        #print(const['type'])
        if const['type'] == None:
           return(True)
        elif const['type'] in '[ dist-x-axis ]':
           dist = const['ref'] - new_point
           status += [dist[0]  < const['tol']]
        elif const['type'] in '[ dist-y-axis ]':
           dist = const['ref'] - new_point
           status += [dist[1]  < const['tol']]
        elif const['type'] in '[ dist-z-axis ]':
           dist = const['ref'] - new_point
           status += [dist[2]  < const['tol']]
           #print(status)

    return(all(status))

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
    #traj = traj.reshape(-1,3)
    bond, angle, dihedral = 0, 0, 0

    for molecule, positions in traj.items():
      for coords in positions:
        #print(ff[molecule])
        bond  += bonded_pot(ff[molecule], coords)
        #print('get here')
        angle += angle_pot(ff[molecule], coords)
        dihedral += dihedral_pot(ff[molecule], coords)
    
    #print(ff['nonbond_params'])
    vdw =  non_bond_interactions(ff, traj)
    print('returned')
    display=True
    if display:
       for term, name in zip([bond, angle, dihedral, vdw],['bonds', 'angle', 'dihedral', 'vdw']):
           print(name, term)

    return(bond + angle + dihedral + vdw)


def metropolis_monte_carlo(ff, name, start, temp, n_repeat, step_length, max_steps, verbose, env_traj, list_of_constraints, sol):
    try:
        traj = {name:env_traj['DOPE'][0]}
        del env_traj['DOPE'][0]
        print(traj)               
    except (KeyError,TypeError):
        traj = {name:[start]}

    count = 0
    prev_E = 0.0
    rejected = 0
    last_point = start
    print('\n++++++++++ STARTING MONTE CARLO MODULE +++++++++\n')
    while count < n_repeat -1:
       print('----------------')     
       while True:
          vector_bundel = norm_sphere()
       #   print(last_point)
          while True:
                new_coord, index = take_step(vector_bundel, step_length, last_point)
                if not is_overlap(new_coord, traj[name], 0.43*0.8,1):
                   break
                else:
                   print('gp')
              

          new_traj = {name:[np.append(traj[name], np.array([new_coord])).reshape(-1,3)]}
          #print(new_traj)
          #exit()

          if len(env_traj) != 0:
             sol_traj_temp = remove_overlap(new_coord, env_traj[sol], 0.47*0.8, sol)                    
             new_traj.update({sol: sol_traj_temp})
             [ new_traj.update({key:values}) for key, values in env_traj.items() if key != sol ]
          
          if constraints(new_coord, list_of_constraints):
             total_E  = Hamiltonion(ff, new_traj, verbose)

             if accaptable(total_E, temp, prev_E):
               if verbose:
                  print('accapted')
                  print(total_E * len(new_traj))
               prev_E = total_E
               traj = new_traj
               last_point = new_coord
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
               print('+++++++++++++++++++++ FATAL ERROR ++++++++++++++++++++++++++++\n')
               print('Exceeded maximum number of steps in the monte-carlo module.')
               print('If you know what you do set -maxconstraints to -1')
          else:
             if verbose:
                print('rejected')         
             rejected = rejected + 1
             vector_bundel = np.delete(vector_bundel, index, axis=0)

    print('++++++++++++++++ RESULTS FORM MONTE CARLO MODULE ++++++++++++++\n')
    print('Total Energy:', Hamiltonion(ff, new_traj, verbose))
    #print('Radius of Gyration:', radius_of_gyr(traj))
    print('Number of rejected MC steps:', rejected)       
    return(traj)

def generate_chains(ff, start, step_length, nexcl, chain, mc_options, env_traj, constraints, sol):
    n_mon, temp, max_steps, verbose, name = mc_options
    traj = start

    if chain in '[ linear-random ]':
       traj = metropolis_monte_carlo(ff, name, start, temp, n_mon, step_length, max_steps, verbose, env_traj, constraints, sol)
    elif chain in '[ PEGylated-bilayer ]':
       traj = metropolis_monte_carlo(ff, name, start, temp, n_mon, step_length, max_steps, verbose, env_traj, constraints, sol)
    elif chain in '[ Polystyrene ]':
       print('PS module not implemented yet!')
       exit()
    return(traj)

def build_system(top_options, env_options, mc_options, outfile):
    topfile, structure_file, chain, conv = top_options 

    conf = read_conf_file(structure_file, 'gro')
    step_length, size, nexcl = 0.322, 0.43, 1   #determine_step_legnth(conf, [0])

    print(env_options)
    if env_options[0] in '[ vac ]':
       env_traj = []
       ff, system = read_top(topfile)
       traj = generate_chains(ff, np.array([0,0,0]), step_length, nexcl, chain, mc_options, env_traj, [{'type':None}], None) 

    elif env_options[0] in '[ bulk ]':
       options = 0
       print('Not active yet!')
       exit()
       env_traj, constraints = create_environment(env_type, options)
       traj = generate_chains(ff, conf, step_length, nexcl, chain, mc_options, env_traj, None)

    elif env_options[0] in '[ bilayer ]':
       env_traj, constraints, head = create_environment(env_options)
       ff, system = read_top(topfile)
       traj = generate_chains(ff, head, step_length, nexcl, chain, mc_options, env_traj, constraints, 'W')
 
    write_gro_file(traj,outfile,ff)
    return(None)   
