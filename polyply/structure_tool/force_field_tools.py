import random
import math
import argparse
import itertools
import numpy as np
import scipy.optimize as opt
import scipy.spatial
import scipy 
from numpy import sqrt, pi, cos, sin, dot, cross, arccos, degrees
from numpy import float as nfl
from numpy.linalg import norm
from string import digits
import multiprocessing
#from polyply.structure_tool.mc_poly_growth import *
from polyply.structure_tool.analysis_funtions import *
from polyply.structure_tool.geometrical_functions import *
from polyply.structure_tool.force_field_tools import *
from multiprocessing import Pool
import time
import scipy.spatial
import polyply.classes.potentials as potentials

global kBa
kb = 1.38964852 * 10**(-23.0) *10**-3.0 # kJ/K


def write_gro_file(top, traj,name,file_name):
    n = len(traj.atom_info)
    out_file = open(file_name, 'w')
    out_file.write('Monte Carlo generated'+name+ '\n')
    out_file.write('{:>3s}{:<8d}{}'.format('',n,'\n'))
    count = 0
    resnum = 1
    print(len(traj.atom_info))
    print(len(traj.positions))
    for pos, info in zip(traj.positions,traj.atom_info):
        #print( top.molecules[name].atoms[info[3]-1].parameters[2])
        res_index = int(top.molecules[name].atoms[info[3]-1].parameters[1])
        res_name = top.molecules[name].atoms[info[3]-1].parameters[2]
        out_file.write('{:>5d}{:<5s}{:>5s}{:5d}{:8.3F}{:8.3F}{:8.3F}{}'.format(res_index,res_name,info[2],info[-1],pos[0],pos[1],pos[2],'\n'))

    out_file.write('{:>2s}{:<.5F} {:<.5F} {:<.5F}'.format('',float(traj.box[0]), float(traj.box[1]), float(traj.box[2])))
    out_file.close()

def get_positions(traj, first_atom, last_atom,centers):
    positions = [ traj.positions[i-1] for i in np.arange(first_atom,last_atom,1) if i < len(traj.positions) ]
    try:
      atom_pos = [ positions[int(i)-1] for i in centers]
      return (atom_pos, True)
    except IndexError:
      return([],False)

def bonded_potential(traj,top,excluded=[None]):
    energies = {}
    atom_count = 0

    for mol, n_mol in top.composition:
       if mol not in excluded:
       
          term_names = [a for a in dir(top.molecules[mol]) if not a.startswith('__') and not callable(getattr(top.molecules[mol],a)) and 
                        a not in ['excl', 'excl_14_list', 'excl_list', 'mol_graph', 'name', 'atoms', 'constraints']]

          for nth_mol in np.arange(0,n_mol,1):
             
              first_atom = atom_count 
              last_atom  = atom_count + len(top.molecules[mol].atoms)
              atom_count = atom_count + len(top.molecules[mol].atoms)
              for term_name in term_names:
                  energy = 0
                  terms = getattr(top.molecules[mol],term_name)
                  
                  for term in terms:
                      positions, proceed = get_positions(traj, first_atom, last_atom,getattr(term,'centers'))                     
                      #print(term.centers, proceed)

                      if proceed:
                         energy += getattr(potentials,getattr(term,'potential'))(term,positions)
                         #print(term.centers)
                
                  try:
                     term_energy = energy + energies[term_name]
                  except KeyError:
                     term_energy = energy
                
                  energies.update({term_name:term_energy})           
       else:
           atom_count = atom_count +  len(top.molecules[mol].atoms) * n_mol    
           #print(atom_count)
    return(energies)
           
def lookup_interaction_parameters(top, atom_A, atom_B, key):

      try:
          coefs = getattr(top,key)[(atom_A, atom_B)].parameters
          pot_form = getattr(top,key)[(atom_A, atom_B)].potential 

      except KeyError:
          coefs = getattr(top,key)[(atom_B, atom_A)].parameters
          pot_form = getattr(top,key)[(atom_B, atom_A)].potential

      coef_A, coef_B = coefs

      if top.defaults['LJ'] != 2:
         if float(coef_A) != 0:
            sigma = (float(coef_B)/float(coef_A))**(1/6)
         else:
            sigma = 0
      else:
         sigma = coef_B  
      return(coefs, sigma, pot_form)

def are_bonded_exception(atom_A, atom_B, molecule, top, key):
    # Note that lists i.e. molecule.bonds are 0 indexed while 
    # the exclusions are indexed starting with 1 as they 
    # are derived from the bond centers which start at 1

    # one particle molecules need to be excepted somehow

   # try:
        #print(top.molecules[molecule].name)
        #print('--',atom_B , getattr(top.molecules[molecule],key)[str(atom_A)] )
        return(any([ atom_B   == atom for atom in getattr(top.molecules[molecule],key)[str(atom_A)] ]))
 #   except KeyError:
  #      return(False)

def nonbonded_potential(dist_matrix, top, softness, eps, verbose):
    print(dist_matrix)
    LJ_energy = 0
    COUL_energy=0
    verbose=False
    if verbose:
       print("using a softness of",softness)

    for item in dist_matrix:
        
        mol_name_A, mol_index_A, atom_type_A, atom_index_A, atom_index_total = item[1]
        mol_name_B, mol_index_B, atom_type_B, atom_index_B, atom_index_total = item[2]
        dist = item[0]
    
        charges = [top.molecules[mol_name_A].atoms[atom_index_A-1].parameters[5],
                   top.molecules[mol_name_B].atoms[atom_index_B-1].parameters[5]]

        coefs, sigma, pot_form = lookup_interaction_parameters(top, atom_type_A, atom_type_B, 'nonbond_params')
        pot_form_COUL = 'Coul'
       
        if mol_name_A == mol_name_B and mol_index_A == mol_index_B:
           if not are_bonded_exception(atom_index_A, atom_index_B, mol_name_A, top,'excl_list'):

              if dist > sigma * softness:
                 LJ_energy   = LJ_energy   + getattr(potentials,pot_form)(coefs, dist,top.defaults['LJ'])
                 COUL_energy = COUL_energy + getattr(potentials,pot_form_COUL)(charges, dist, eps)
 #                print(item)
              else:
                 if verbose:
                    print(mol_name_A)
                    print('overalp - category I')
                 return(math.inf, math.inf)

           elif are_bonded_exception(atom_index_A, atom_index_B, mol_name_A, top,'excl_list'):
                LJ_energy = LJ_energy + 0
                COUL_energy = COUL_energy + 0

           elif are_bonded_exception(atom_index_A, atom_index_B, mol_name_A, top,'excl_14_list'):
                try:
                   coefs, sigma, pot_form = lookup_interaction_parameters(top, atom_type_A, atom_type_B, 'nonbond_14_pairs')
                except AttributeError:
                   if verbose:
                      print('No 1-4 interactions found')

                if dist > sigma * softness: 
                   LJ_energy = LJ_energy + getattr(potentials,pot_form)(coefs, dist,top.defaults['LJ'])
                   COUL_energy = COUL_energy + getattr(potentials,pot_form_COUL)(charges, dist, eps)
                else:
                   if verbose:
                      print('A-self-overlap-B')
                      print(dist)
                   return(math.inf, math.inf)
        else: 

           if dist > sigma * softness:
              LJ_energy   = LJ_energy   + getattr(potentials,pot_form)(coefs, dist, top.defaults['LJ'])
              COUL_energy = COUL_energy + getattr(potentials,pot_form_COUL)(charges, dist, eps)
              #print(item)
           else:
              if verbose:
                 print(item)
                 print('A-self-overlap-C')
                 print(dist)
              return(math.inf, math.inf)  
        
    return(LJ_energy, COUL_energy)
