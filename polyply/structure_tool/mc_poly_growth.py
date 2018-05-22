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


def find_central_starting_point(coordinates):
    x_coords = [ atom[0] for atom in coordinates ]
    y_coords = [ atom[1] for atom in coordinates ]
    z_coords = [ atom[2] for atom in coordinates ]
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
    display = True   
    bonded_energies= bonded_potential(traj, top, [sol])
    vdw, coulomb = nonbonded_potential(traj.dist_matrix, top, softness, eps, display)
    print(bonded_energies)
    if display:
       for item in bonded_energies:
           print(item)
       print('vdw:',vdw)
       print('electrostatic',coulomb)
   
    bond = sum([ i[0] for i in bonded_energies  ])
    return( bond + vdw + coulomb)

def minimize(new_traj, name):
    pass

##############################################################################################################################
#
#                                             The Pseudopolis Monte-Carlo Method
#
##############################################################################################################################

def constraint_vectors(constraints):
    vectors = norm_sphere()
    return(vectors)


def look_up_reference_pos(traj, top, n_atoms_placed, name):

    # based on the graph spanned by the bonds we look up all directly bonded atoms
    attached_atoms_indices = top.molecules[name].neighborhood(n_atoms_placed+1, 1)
  #  print('molecule:',[pos for pos, mol in zip(traj.positions, traj.atom_info) if mol[0] == name])
 #   print('placing atom:',n_atoms_placed+1)
 #   print('attached atoms',attached_atoms_indices)  
    # next we select all those that already have been placed, note that the traj is zero indexed while the bond-list is not
    ref_atom_index = [ atom for atom in attached_atoms_indices if atom <= n_atoms_placed ][0] 
 #   print(ref_atom_index) 
 #   print([pos for pos, mol in zip(traj.positions, traj.atom_info) if mol[0] == name and mol[3] == ref_atom_index])
    # finally we select the reference position as theatom which is part of the last molecule in the trajectory
    ref_position = [  pos for pos, mol in zip(traj.positions, traj.atom_info) if mol[0] == name and mol[3] == ref_atom_index  ][-1]
    return ref_position

def determine_step_length(name, top, traj):
    # first we look how many atoms of the same moleculetype exist
    n_atoms_mol = len([ info[0]  for info in  traj.atom_info if info[0] == name])
    
    # then we check how many complete molecules we have and substract them from the total
    n_atoms_placed = n_atoms_mol - (n_atoms_mol // len(top.molecules[name].atoms)) * len(top.molecules[name].atoms)

    # based on the number of atoms of the incomplete molecule we select the reference position
    ref_position = look_up_reference_pos(traj, top, n_atoms_placed, name)

    # the current atom is the number of placed atoms incremented by 1; note we place all atoms consecutively
    step_length =  float(top.molecules[name].bonds[n_atoms_placed+1].parameters[0])
    print('n_atoms_placed:',n_atoms_placed)
    return(step_length, ref_position, n_atoms_placed+1)


def attempt(traj, top, step_length, ref_position, vectors, name, mol_index, mol_atom_index, sol_name, softness, eps, cut_off):

    # a new point is generated by randomly selecting a vector
    index     = random.randint(0, len(vectors) - 1)
    new_point = ref_position + vectors[index] * step_length
    
    new_traj             = trajectory('temp')
    new_traj.positions   = traj.positions
    new_traj.atom_info[:]   = traj.atom_info
    new_traj.box         = traj.box

    # the new point is added to the trajectory and all overlapping molecules are removed
    new_traj.add_atom(name, mol_index, mol_atom_index,new_point,len(traj.positions), top)
    new_traj.distance_matrix(cut_off,top)
    new_traj.remove_overlap(top, name, sol_name, cut_off, softness)

    # the energy of this cleaned trajectory is computed
    new_energy = Hamiltonion(top, new_traj, False, softness, eps,sol_name)

    return(new_energy, new_traj)
           
def take_pseudo_step(traj, top, maxsteps, name, sol_name, verbose, softness, eps,sol, mol_index, n_min_steps, cut_off,step_count,prev_energy,temp,rejected,max_steps):

    step_length, ref_position, mol_atom_index = determine_step_length(name, top, traj)
    vectors = constraint_vectors(constraints)    
    subcount = 0

    while True:

         print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
         new_energy, new_traj = attempt(traj, top, step_length, ref_position, vectors, name, mol_index, mol_atom_index, sol_name, softness, eps, cut_off)
      

 #        print(new_traj.positions)
#         print(traj.positions)
         print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

         if step_count//n_min_steps - step_count == 0:
            return(new_energy, new_traj, rejected, True)

         elif accaptable(new_energy,temp, prev_energy):
            return(new_energy, new_traj, rejected, False)

         elif subcount < max_steps:
              if verbose:
                 print(new_energy)
                 print('rejected')
              
              rejected = rejected + 1
              subcount = subcount + 1
              vectors = vectors[subcount:]
         else:
              print('+++++++++++++++++++++ FATAL ERROR ++++++++++++++++++++++++++++\n')
              print('Exceeded maximum number of steps in the monte-carlo module.')
              print('If you know what you do set -maxsteps to -1')


def pseudopolis_monte_carlo(top, traj, start, name,  temp, max_steps, verbose, list_of_constraints, sol, cut_off, eps, softness,n_min_steps):
 
###### 1. Check if it is a restart 
 
  
    n_repeat = len(top.molecules[name].atoms)
    n_mol = sum([entry[1] for entry in  top.composition if entry[0] == name])-1

    # yes it is a restart; we need to look-up how many atoms already were placed
    if start == None:
       print("Restarting from given structure!")
       n_atoms_mol = len([ info[0]  for info in  traj.atom_info if info[0] == name])
       count = n_atoms_mol - (n_atoms_mol / len(top.molecules[name].atoms)) * len(top.molecules[name].atoms) 
       
    # no it is a new start; but we have solvent so we start from the center of the box
    elif start == 'center':
       print("Will start from center of box!")
       ref_point = find_central_starting_point(traj.positions) 
       traj.add_atom(name, n_mol, 1, ref_point, len(traj.positions), top)
       traj.remove_overlap(top, name, sol, cut_off, softness)
       count = 1

    else:
       print("Will start from scratch! Hold tight!")
       traj.add_atom(name, n_mol, 0,np.array([0,0,0]),len(traj.positions), top)
       count = 1

###### 2. Construct the inital distance matrix of the entire box       
  
    print(len(traj.positions))
    traj.distance_matrix(cut_off,top)
    print(len(traj.positions))
    if verbose:
       print('computed ',len(traj.dist_matrix),'distance pairs.')


###### 3. Compute the intial energy

    prev_energy = Hamiltonion(top, traj, verbose, softness, eps,sol)
    
###### 4. loop over n_repeats
    n_min_steps = n_repeat + 1
    rejected=0
    print('\n+++++++++++++++ STARTING PSEUDOPOLIS-MONTE CARLO MODULE ++++++++++++++\n')
    while count < n_repeat:
 
          energy, new_traj, new_rejected, minimize = take_pseudo_step(traj, top, max_steps, name, sol,verbose, softness, eps,sol,n_mol,n_min_steps,cut_off,count,prev_energy,temp,rejected,max_steps)   
          print('~~~~~~~~~~~~~',count,'~~~~~~~~~~~~~~~~')
          if minimize:
             traj, energy = minimize(new_traj, name)
          else:
             traj = new_traj
            
          rejected = rejected + new_rejected
          count = count + 1

    print('++++++++++++++++ RESULTS FORM MONTE CARLO MODULE ++++++++++++++\n')
    print('Total Energy:', Hamiltonion(top, new_traj, verbose, softness, eps,sol))
    print('Number of rejected MC steps:', rejected)       
    return(traj)


################################################################################################################
#                                 you're leaving the pseudopolis MC 
################################################################################################################

def build_system(top_options, env_options, mc_options, outfile, magic_numbers):
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
    start='center'
    n_min_steps=0
    traj = pseudopolis_monte_carlo(top, traj, start, name,  temp, max_steps, verbose, [], sol, cut_off, eps, softness, n_min_steps)
  
    write_gro_file(top,traj,name,'test.gro')
    return(None)   
