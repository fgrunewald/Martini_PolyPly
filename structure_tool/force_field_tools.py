import random
import math
import argparse
import itertools
import numpy as np
import scipy.optimize as opt
from numpy import sqrt, pi, cos, sin, dot, cross, arccos, degrees
from numpy import float as nfl
from numpy.linalg import norm
from string import digits
import multiprocessing
#from Martini_PolyPly.structure_tool.mc_poly_growth import *
from Martini_PolyPly.structure_tool.analysis_funtions import *
from Martini_PolyPly.structure_tool.geometrical_functions import *
from Martini_PolyPly.structure_tool.force_field_tools import *
from multiprocessing import Pool
import time

global kBa
kb = 1.38964852 * 10**(-23.0) *10**-3.0 # kJ/K

def read_top(name):
    itpfiles =  []
    system, ff = {}, {}
    empty = 0
    section = 'Random'
    with open(name) as f:
         lines = f.readlines()
         for line in lines:
           #print(line)
           if len(line.replace('\n', '').split()) == 0:
              empty = empty +  1
           elif not any([ word in ';' for word in line.split()]):
              if any([ word in '[ [ ]' for word in line.split()]):
                 section = line.replace('\n', '').split()[1]
                 print(section)
              elif section in '[ molecules ]':
                 name, n_mol = line.replace('\n', '').split()
                 system.update({name: int(n_mol)})
              elif section in '[ System ]':
                 print('Reading in', line)  
              elif any([ word in '#include' for word in line.split()]):
                 itpfiles += [line.replace('\n', ' ').replace('\"', ' ').split()[1]]
                 print(line)

    for itp in itpfiles:
        name, parameters =  read_itp(itp)
        ff.update({name: parameters})
       
    return(ff, system)

def read_itp(name):
    print('Reading',name)
    ff = {}
    molecules, atoms, bonds, angles, dih, constraints, virtual_sitsn = [], [], [], [], [], [], []
    nonbond_params = {}
    section = 'moleculetype'
    empty=0
    with open(name) as f:
         lines=f.readlines()
         for line in lines:
           if len(line.replace('\n', '').split()) == 0:
              empty = empty +  1
           elif not any([ word in ';' for word in line.split()]):
 
              if any([ word in '[ [ ]' for word in line.split()]):
                 section = line.replace('\n', '').split()[1]
 
              elif section in '[ moleculetype ]':
 #                  print(line)
                   name, nexcl = line.replace('\n', '').split()
                   molecules.append({'name':name,'nexcl':nfl(nexcl)})
                   molecule_type = name
                    
              elif section in '[ atoms ]':
                  n, typ, resnr, res, atom, cgnr, charge, mass = line.replace('\n', '').split()
                  atoms.append({'n': int(n), 'typ':typ ,'atom':atom, 'charge':nfl(charge)})
 
              elif section in '[ nonbond_params ]':
                  atom1, atom2, f, sigma, epsilon = line.replace('\n', '').split()
                  nonbond_params.update({(atom1, atom2): {'sigma':nfl(sigma), 'epsilon':nfl(epsilon)}})
 
              elif section in '[ bonds ]':
                  A, B, f, ref, k0 = line.replace('\n', '').split()
                  bonds.append({'pairs':[int(A),int(B)], 'k0':nfl(k0), 'ref':nfl(ref)})
 
              elif section in '[ angles ]':
                  A, B, C, f, ref, k0 = line.replace('\n', '').split()
                  angles.append({'pairs': [int(A), int(B), int(C)], 'k0':nfl(k0), 'ref':nfl(ref)})
  
              elif section in '[ dihedrals ]':
                  A, B, C, D, f, ref, k0, n = line.replace('\n', '').split()
                  dih.append({'pairs':[int(A),int(B),int(C), int(D)], 'k0':nfl(k0), 'f':nfl(f), 'n':nfl(n), 'ref':nfl(ref)})
  
              elif section in '[ constraints ]':
                  A, B, f, ref = line.replace('\n', '').split()
                  constraints.append({'pairs':[A, B], 'f':nfl(f), 'ref':nfl(ref)})
  
    if len(nonbond_params) != 0:
       return('nonbond_params', nonbond_params)
    else:
       return(molecule_type,{'nexcl':nfl(nexcl), 'atoms':atoms, 'bonds':bonds, 'angles':angles, 'constraints':constraints, 'dih':dih})

def convert_constraints(ff):
    print('++++++++++++++++++++++ Converting Constraints +++++++++++++++++++++++')
    for molecule in ff:
     if molecule != 'nonbond_params':
      if len(ff[molecule]['constraints']) != 0:
        new_bonds = []
        #print(len(ff[molecule]['bonds']))
        new_bonds = [ {'pairs':[int(term['pairs'][0]), int(term['pairs'][1])], 'ref':term['ref'], 'k0': nfl(9000.0)} for term in ff[molecule]['constraints'] ]
        new_bonds =  ff[molecule]['bonds'] + new_bonds
        ff[molecule].update({'bonds':new_bonds})
    return(ff)

def write_gro_file(data, name, ff, box):
    n = sum([ len(coords) for resname, list_of_coords in data.items() for coords in list_of_coords ])
    out_file = open(name, 'w')
    out_file.write('Monte Carlo generated PEO'+'\n')
    out_file.write('{:>3s}{:<8d}{}'.format('',n,'\n'))
    count = 0
    resnum = 1
  
    for resname, mols in data.items():
       for coords in mols:
         resnum = resnum + 1
         for index, line in enumerate(coords):
             atomtype = ff[resname]['atoms'][index]['atom']
             count = count + 1
             out_file.write('{:>5d}{:<5s}{:>5s}{:5d}{:8.3F}{:8.3F}{:8.3F}{}'.format(resnum, resname, atomtype, count, line[0], line[1], line[2],'\n'))
    out_file.write('{:>2s}{:<.5F} {:<.5F} {:<.5F}'.format('',float(box[0]), float(box[1]), float(box[2])))
    out_file.close()
    return(None)

def write_log_file(infos, fomat_info):
    log_file = open('logfile.log', 'w')
    for line in infos:
        log_file.write(format_info.format(info))
    log_file.close()
    return(None)




def pot_I( val ,k0, ref):
    return(0.5 * k0 * (val-ref)**2.0)

def LJ(C6, C12, r):
  #  print(sig)
 #   print(eps)
    #return(4.0 * eps * ( (sig/r)**12.0 - (sig/r)**6.0))
    return( (C12/r**12.0 - C6/r**6.0) )    

def proper_dih(dih, k0, ref, n):
    return( k0 * (1 + cos(n * np.radians(dih) - np.radians(ref))))

def legal(term, traj):
    status_A = all( [index <= len(traj) for index in term['pairs']])
    if status_A:
       coords = [traj[i - 1] for i in term['pairs']]
       #print(coords)
       status_B = all([any(coord != np.array([0,0,0])) for coord in coords])
       return(status_B)
    else:
       return(status_A)

def bonded_pot(ff, traj):
    bond_pairs = [ norm(traj[(term['pairs'][0] - 1)] - traj[(term['pairs'][1] - 1)]) for term in ff['bonds'] if legal(term, traj)]
    return(sum([pot_I(dist, term['k0'], term['ref']) for term, dist in zip(ff['bonds'], bond_pairs)]))

def angle_pot(ff, traj):
    angles =  [ angle(traj[(term['pairs'][0] - 1)], traj[(term['pairs'][1] - 1)], traj[(term['pairs'][2] - 1)]) for term in ff['angles'] if legal(term, traj) ]
    return(sum([pot_I(np.radians(ang),term['k0'], np.radians(term['ref'])) for term, ang in zip(ff['angles'], angles)]))

def dihedral_pot(ff, traj):
    dih_ang = [dih(traj[(t['pairs'][0] - 1)], traj[(t['pairs'][1] - 1)], traj[(t['pairs'][2] - 1)], traj[(t['pairs'][3] -1)]) for t in ff['dih'] if legal(t, traj)]
    return(sum([proper_dih(ang, term['k0'], term['ref'], term['n']) for term, ang in zip(ff['dih'], dih_ang)]))

def are_bonded(atom_A,atom_B,molecule_A,ff):
    atom_A = atom_A + 1
    atom_B = atom_B + 1
    if len(ff[molecule_A]['bonds']) > 1:
      for bond in ff[molecule_A]['bonds']:
        index_A, index_B = bond['pairs'][0], bond['pairs'][1]
        #print(molecule_A)
        #print(index_A, index_B)
        #print(atom_A, atom_B)
        if index_A == atom_A and index_B == atom_B:
           #print(molecule_A)
           #print( index_A == atom_A )
           #print(  index_B == atom_B)
           #exit()
           return(True)
        elif index_A == atom_B and index_B == atom_A:
           return(True)
      else:
           return(False)
    else:
      return(False)


def partition(lst, n):
    parts = len(lst) / n
    return([lst[round(parts * i):round(parts * (i + 1))] for i in range(n)])

def non_bond_interactions(ff, traj):
    start = time.time()
    #print('------> computing non-bonded interactions')
    '''
    The outermost loop over molecules is computed in parallel. To do so we, however,
    have to flaten the traj, which for not so long trajs is OK, as long as they can
    fit in the RAM.
    '''
    flat_traj = [ [ molecule, pos ]  for molecule, pos_mols in traj.items() for pos in pos_mols ]
    partitioned_data = partition(flat_traj, 4)
    data = [ (ff, part, traj) for part in partitioned_data ]
    with  multiprocessing.Pool(4) as p:
          energy = p.map(Vdw_pot, data)
    energy = sum(energy)
    return(energy)

def Vdw_pot(input_data):
    ff, traj_part, traj = input_data
    energy=0
    trying=0
    for coords_A in traj_part:   
          molecule_A = coords_A[0]
          pos_mol_A = coords_A[1]
       #   print(coords_A)
          for i, point_A in enumerate(pos_mol_A):
             for molecule_B, coords_B in traj.items():
               for pos_mol_B in coords_B:
                  for j, point_B in enumerate(pos_mol_B):
                     start = time.time()
                     dist = norm(point_A - point_B)
                     trying = trying + (time.time()-start)
                     atom_A, atom_B = ff[molecule_A]['atoms'][i]['typ'], ff[molecule_B]['atoms'][j]['typ']
                     try:
                         epsilon = ff['nonbond_params'][(atom_A, atom_B)]['epsilon']
                         sigma = ff['nonbond_params'][(atom_A, atom_B)]['sigma']
                     except KeyError:
                         epsilon = ff['nonbond_params'][(atom_B, atom_A)]['epsilon']
                         sigma = ff['nonbond_params'][(atom_B, atom_A)]['sigma']
                     if molecule_A == molecule_B:
                        if i-j != 0:
                          if not are_bonded(i, j, molecule_A, ff): #> ff[molecule_A]['nexcl']:
                           energy = energy + LJ(sigma, epsilon, dist)                           
                           if dist < 0.75*0.43:
                              energy = math.inf
                            #  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
                              print('Self-Overlap')
                              print(dist)
                            #  print('atomA:',atom_A)
                            #  print('atomB:',atom_B)
                            #  print(i,j)
                            #  print('diff:',j-i)
                            #  print(are_bonded(i, j, molecule_A, ff))
                            #  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
                              return(energy)
                        else:
                           energy = energy   
                     else:
                        energy = energy + LJ(sigma, epsilon, dist)
                        #if dist < 0.5:
                           #print(LJ(sigma, epsilon, dist))
                        if dist < 0.75*0.47:
                           energy = math.inf
                         #  print('Overlap')
                         #  print(molecule_B)
                         #  print(molecule_A)
                           return(energy)
   # print('Spengt',trying,'seconds')
    return(energy)

