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
from polyply.structure_tool.mc_poly_growth import *
from polyply.structure_tool.analysis_funtions import *
from polyply.structure_tool.geometrical_functions import *
from polyply.structure_tool.force_field_tools import *
from polyply.structure_tool.environment import *
from polyply.classes.topology_class import *
from polyply.classes.traj_class import *
from polyply.classes.topology_class import *
from polyply.classes.traj_class import *
from polyply.classes.potentials import *

#######################################################################################################
#
#                       Functions for the pseudo Monte-Carlo Method
#
#######################################################################################################


def find_central_starting_point(coordinates, sol):
    x_coords = [ atom[0] for sol_name in sol for mol in coordinates[sol_name] for atom in mol ]
    y_coords = [ atom[1] for sol_name in sol for mol in coordinates[sol_name] for atom in mol ]
    z_coords = [ atom[2] for sol_name in sol for mol in coordinates[sol_name] for atom in mol ]
    point_A = np.array([max(x_coords), max(y_coords), max(z_coords)])
    point_B = np.array([min(x_coords), min(y_coords), min(z_coords)])
    diagonal = point_A - point_B
    starting_point = point_B + 0.5 * diagonal
    return(starting_point)

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

def is_overlap(new_point, traj, tol, bonds, current):
    bonds = [ index - 1 for pair in bonds for index in pair['pairs'] if index != current + 1 ]
    for index, coords in enumerate(traj):
      if index not in bonds:
        for point in coords:
           if norm(point - new_point) < tol:
              return(True)
    return(False)


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
    #print(status)

    return(all(status))

def Hamiltonion(top, traj, display, softness, eps,sol):
    
    bonded_energies= bonded_potential(traj, top, [sol])
    vdw, coulomb = nonbonded_potential(traj.dist_matrix, top, softness, eps, display)
   
    if display:
       for item in bonded_energies:
           print(item)
       print('vdw':,vwd)
       print('electrostatic',coulomb)

    return(vdw + coulomb)



##############################################################################################################################
#
#                                             The Pseudopolis Monte-Carlo Method
#
##############################################################################################################################


def attempt(traj, top, step_length, vectors):

    last_point = traj.positions[-1]
    index = random.randint(0, len(vectors) - 1)
    new_point = last_point + vectors[index] * step_length
    traj.positions.append(new_point)
    traj.remove_overlaps(name)
    new_energy = Hamiltonion(traj,top)

    return(new_energy, traj)
           
def take_step():

    step_length = determine_step_length()
    vectors = constraint_vectors()    

    while True:

         new_energy, new_traj = attempt(traj, top, step_length, vectors)
   
         if type(step_count/min_steps) == int:
            return(new_energy, new_traj)

         elif accaptable(new_energy, prev_energy):
            return(new_energy, new_traj)

         elif subcount < max_steps:
              if verbose:
                 print(total_E)
                 print('rejected')

               vector_bundel = vector_bundel[]
               rejected = rejected + 1
               subcount = subcount + 1


   if len(env_traj) != 0:
       sol_traj_temp = remove_overlap(new_coord, env_traj[sol], 0.43*0.80, sol)                    
       new_traj.update({sol: sol_traj_temp})
       [ new_traj.update({key:values}) for key, values in env_traj.items() if key != sol ]
 
       dist_mat = construct_dist_mat(new_traj, cut_off, start=True)



    return(new_item, index)


def pseudo_step():

    while True:
          vector_bundel = norm_sphere()
          new_coord, index, new_E = take_step(vector_bundel, step_length, last_point)

          if accaptable(new_E, temp, prev_E):
             if verbose:
                print('accapted')
             prev_E = total_E
             last_point = new_coord
             count = count + 1
             break

          elif subcount < max_steps:
               if verbose:
                  print(total_E)
                  print('rejected')         
               
               vector_bundel = vector_bundel[]
               rejected = rejected + 1
               subcount = subcount + 1

             else:
               print('+++++++++++++++++++++ FATAL ERROR ++++++++++++++++++++++++++++\n')
               print('Exceeded maximum number of steps in the monte-carlo module.')
               print('If you know what you do set -maxsteps to -1')
               return(traj)
          else:
             if verbose:
                print('rejected')         
             rejected = rejected + 1
                   vector_bundel = np.delete(vector_bundel, index, axis=0)



def pseudopolis_monte_carlo(top, name, traj, start,temp, max_steps, verbose, list_of_constraints, sol, cut_off, eps, softness):
 
###### 1. Check if it is a restart 
         #    -> set the number of iterations n_repeat
         #    -> set the last point 

    if start != None:
       print("Restarting from given structure!")
       last_point = traj.mol_pos[name][-1]
       n_repeat = len(top.molecules[name].atoms) - len(traj.mol_pos[name])
    elif start == 'center':
       print("Will start from center of box!")
       last_point = find_central_starting_point(traj.mol_pos, sol) 
       n_repeat = len(top.molecules[name].atoms)
    else:
       print("Will start from scratch! Hold tight!")
       last_point = np.array([0,0,0])
       n_repeat = len(top.molecules[name].atoms)

###### 2. Construct the inital distance matrix of the entire box       
  
    traj.distance_matrix(cut_off)
    if verbose:
       print('computed ',len(traj.dist_matrix),'distance pairs.')


###### 3. Compute the intial energy

    prev_E = Hamiltonion(top, traj, verbose, softness, eps,sol)
    
###### 4. loop over n_repeats

    rejected=0
     
    print('\n+++++++++++++++ STARTING PSEUDOPOLIS-MONTE CARLO MODULE ++++++++++++++\n')
    while count < n_repeat:
 
    
    print('~~~~~~~~~~~~~~~~',count,'~~~~~~~~~~~~~~')
    step_length, ref_coord, bonded = determine_step_length(ff, count, traj[name][0], name, start, offset)
    last_point = ref_coord
    subcount=0     

    print('++++++++++++++++ RESULTS FORM MONTE CARLO MODULE ++++++++++++++\n')
    print('Total Energy:', Hamiltonion(ff, traj, dist_mat, True, eps, cut_off, form, softness))
    #print('Radius of Gyration:', radius_of_gyr(traj))
    print('Number of rejected MC steps:', rejected)       
    return(traj)


################################################################################################################
#                                 you're leaving the pseudopolis MC 
################################################################################################################

def build_system(top_options, env_options, mc_options, outfile, magic_numbers):
    # some magic numbers which should go to input level at some point
    cut_off, softness, eps,  verbose = magic_numbers
    topfile = top_options
    temp, max_steps, verbose, name = mc_options
    sol, sysfile, start = env_options 

    top_format = topology_format('gromacs','topology_format.txt')
    top = topology.from_gromacs_topfile(topfile,top_format)

    if sysfile != None:
       traj = trajectory.from_gro_file(sysfile,top)
    else:
       traj = trajectory(name)  
       traj.add_molecule(name, positions=np.array([0,0,0]), top=top)
   
    traj = pseudopolis_monte_carlo(top, name, traj, start, temp, max_steps, verbose, [{'type':None}],  None,  cut_off, eps, softness)
  
    write_gro_file(traj,outfile,ff, box)
    return(None)   
