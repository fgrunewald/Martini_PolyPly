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
 #   print(new_point)  
 #   print(list_of_constraints)
    for const in list_of_constraints:
        #print(const['type'])
        if const['type'] == None:
           return(True)
        elif const['type'] in '[ dist-x-axis ]':
           dist = const['ref'] - new_point
           status += [abs(dist[0])  < const['tol']]
        elif const['type'] in '[ dist-y-axis ]':
           dist = const['ref'] - new_point
           status += [abs(dist[1])  < const['tol']]
        elif const['type'] in '[ dist-z-axis ]':
           dist = const['ref'] - new_point
           status += [dist[2]  < const['tol']]
    print(status)

    return(all(status))

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
    #print('returned')
    display=False
    if display:
       for term, name in zip([bond, angle, dihedral, vdw],['bonds', 'angle', 'dihedral', 'vdw']):
           print(name, term)

    return(bond + angle + dihedral + vdw)

def is_in_pair_a(pair, value): 
    #print(pair)
    if value == pair[0] and value > pair[1]:
       return(True)
    elif value == pair[1] and value > pair[0]:
       return(True)    
    else:
       return(False)

def is_in_pair_b(pair, value):
    #print(pair) 
    #print(value)
    if value == pair[0]:
       return(True)
    elif value == pair[1]:
       return(True)    
    else:
       return(False)

def is_in_pair_c(pair, n):
    print(pair)
    for i in np.arange(0,n,1):
        if i == pair[0]:
           return(True)
        elif i == pair[1]:
           return(True)    
    else:
       return(False)

def determine_step_length(ff, count, traj, name, start, offset):
    count = count + offset
    current_atom = ff[name]['atoms'][count]['n']
    prev_atom = ff[name]['atoms'][count - 1 ]['n']
  
    bonds_a = [ bond  for bond in ff[name]['bonds'] if is_in_pair_a(bond['pairs'], current_atom)]
    bonds_b = [bond for bond in bonds_a if is_in_pair_b(bond['pairs'], prev_atom)]
   
    if len(bonds_b) == 0:
          bonds_c = [ bond for bond in bonds_a if is_in_pair_c(bond['pairs'], len(traj))]
          index = [n for n in bonds_c[0]['pairs'] if n != count ][0] -1 
          ref_coord = traj[index]
          step_length = bonds_c[0]['ref']
    else:
          step_length = bonds_b[0]['ref']
          ref_coord = traj[count - 1 ]
    return(step_length, ref_coord)


def metropolis_monte_carlo(ff, name, start, temp, n_repeat, max_steps, verbose, env_traj, list_of_constraints, sol, offset):
 
    try:
        traj = {name:[env_traj['DOPE'][0]]}
        del env_traj['DOPE'][0]
        #print(traj)               
    except (KeyError,TypeError):
        traj = {name:[start]}
 
    count = 0
    print('\n++++++++++ STARTING MONTE CARLO MODULE +++++++++\n')
    while count < n_repeat:
       print('~~~~~~~~~~~~~~~~',count,'~~~~~~~~~~~~~~')
       if count == 0:
          prev_E = 0
          rejected = 0
          count = count + 1
       else:
          print(list_of_constraints)
          step_length, ref_coord = determine_step_length(ff, count, traj[name][0], name, start, offset)
          old = ff[name]['bonds'][count - 1]
          last_point = ref_coord     
          while True:
                vector_bundel = norm_sphere()
                while True:
                      new_coord, index = take_step(vector_bundel, step_length, last_point)
                      if not is_overlap(new_coord, traj[name], 0.43*0.8,1):
                         break
                      #else:
                         #print('gp')
               # print(traj)
                new_traj = {name:[np.append(traj[name],np.array([new_coord])).reshape(-1,3)]}
                  

                if len(env_traj) != 0:
                   sol_traj_temp = remove_overlap(new_coord, env_traj[sol], 0.47*0.8, sol)                    
                   new_traj.update({sol: sol_traj_temp})
                   [ new_traj.update({key:values}) for key, values in env_traj.items() if key != sol ]
                
 
                if constraints(new_coord, list_of_constraints):
                   total_E  = Hamiltonion(ff, new_traj, verbose)
              #     print(new_traj)
                   if accaptable(total_E, temp, prev_E):
                     if verbose:
                        print('accapted')
                        print(total_E * len(new_traj))
                     prev_E = total_E
                     traj = new_traj
               #      print(traj)
                     last_point = new_coord
                  #   print(count)
                     count = count + 1
                     break
                   elif count < max_steps:
                     #print('rejected')
                     #exit()
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

def generate_chains(ff, start, chain, mc_options, env_traj, constraints, sol):
    n_mon, temp, max_steps, verbose, name = mc_options
    traj = start

    if chain in '[ linear-random ]':
       traj = metropolis_monte_carlo(ff, name, start, temp, n_mon, max_steps, verbose, env_traj, constraints, sol, 0)
    elif chain in '[ PEGylated-bilayer ]':
       traj = metropolis_monte_carlo(ff, name, start, temp, n_mon, max_steps, verbose, env_traj, constraints, sol, 11)
    elif chain in '[ Polystyrene ]':
       print('PS module not implemented yet!')
       exit()
    return(traj)

def build_system(top_options, env_options, mc_options, outfile):
    topfile = top_options
 
    size, nexcl = 0.43, 1
    box = np.array([10,10,10])
    #print(env_options)
    if env_options[0] in '[ vac ]':
       env_traj = []
       ff, system = read_top(topfile)
       traj = generate_chains(ff, np.array([0,0,0]), chain, mc_options, env_traj, [{'type':None}], None)

    elif env_options[0] in '[ bulk ]':
       options = 0
       print('Not active yet!')
       exit()
       env_traj, constraints = create_environment(env_type, options)
       traj = generate_chains(ff, conf, step_length, nexcl, chain, mc_options, env_traj, None)

    elif env_options[0] in '[ bilayer ]':
       env_traj, constraints, head, box = create_environment(env_options)
       print(head)
       ff, system = read_top(topfile)
       traj = generate_chains(ff, head, env_options[0], mc_options, env_traj, constraints, 'W')
 
    write_gro_file(traj,outfile,ff, box)
    return(None)   
