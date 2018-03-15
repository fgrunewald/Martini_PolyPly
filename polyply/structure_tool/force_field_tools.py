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
                 #print(section)
              elif section in '[ molecules ]':
                 name, n_mol = line.replace('\n', '').split()
                 system.update({name: int(n_mol)})
              elif section in '[ System ]':
                 q=1
                 #print('Reading in', line)  
              elif any([ word in '#include' for word in line.split()]):
                 itpfiles += [line.replace('\n', ' ').replace('\"', ' ').split()[1]]
                 #print(line)

    for itp in itpfiles:
        name, parameters =  read_itp(itp)
        ff.update({name: parameters})
    print(ff)      
    return(ff, system)

def read_itp(name):
    print('Reading',name)
    ff = {}
    molecules, atoms, bonds, angles, dih, constraints, virtual_sitsn = [], [], [], [], [], [], []
    nonbond_params = {}
    section = 'random'
    empty=0
    with open(name) as f:
         lines=f.readlines()
         for line in lines:
           if len(line.replace('\n', '').split()) == 0:
              empty = empty +  1
           elif not any([ line[i] in ';' for i in np.arange(0,len(line))]):
 
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
              elif section in '[ defaults ]':
                  LJ = int(line.replace('\n', '').split()[0])
                  if LJ == 1:
                     form='C6C12'
                  else:
                     form='sigeps'
                  nonbond_params.update({'functype': form})
 
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

def pot_I( val ,k0, ref):
    return(0.5 * k0 * (val-ref)**2.0)   

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

def construct_bonded_exclusions(ff, nexcl):
    bonded_lists = {}
    for molecule in ff.items() if molecule != 'nonbonded_params':
        mol_graph = construct_mol_graph(molecule['bonds'])
        bonded_list = []
        for atom in molecule['atoms']:
            bonded_list += [(atom['n'], neighborhood(mol_graph, atom['n'], nexcl)]
        bonded_lists.update({molecule:bonded_list})
    return(bonded_lists)

def construct_mol_graph(bonds):
    G = nx.Graph()
    edges = [ (entry['pairs'][0], entry['pairs'][1]) ]
    G.add_edges_from(edges)
    return(G)

def neighborhood(G, node, n):
    '''
     Adobted from: https://stackoverflow.com/questions/22742754/finding-the-n-degree-neighborhood-of-a-node
    '''
    path_lengths = nx.single_source_dijkstra_path_length(G, node)
    neighbours=[ node for node, length in path_lengths.items() if length <= n]
    return(neighbours)


def coulomb(ca, cb, dist,eps):
    return(1/(4*np.pi*eps)*(ca*cb)/dist**2.0)

def LJ(A, B, r, form):
    #print('go here')
    if form == 'C6C12':
       #print(A/r**12.0 - B/r**6.0)
       return( A/r**12.0 - B/r**6.0 )
    elif form == 'sigeps':
       return(4.0 * B * ( (A/r)**12.0 - (A/r)**6.0) )

def nonbonded_potential(dist_mat, ff, softness, eps, form, verbose):
    energy = 0
    e_pot=0
    verbose=True
   # print(len(dist_mat))
   
    for key, dist in dist_mat.items():
        #print(key)
        atom_A, atom_B = ff[key[0]]['atoms'][key[2]]['typ'], ff[key[3]]['atoms'][key[5]]['typ']
        charge_A, charge_B = ff[key[0]]['atoms'][key[2]]['charge'], ff[key[3]]['atoms'][key[5]]['charge']

        try:
            coef_A = ff['nonbond_params'][(atom_A, atom_B)]['epsilon']
            coef_B = ff['nonbond_params'][(atom_A, atom_B)]['sigma']
        except KeyError:
            coef_A = ff['nonbond_params'][(atom_B, atom_A)]['epsilon']
            coef_B = ff['nonbond_params'][(atom_B, atom_A)]['sigma']

        if form == 'C6C12':
           sigma = (coef_A/coef_B)**(1/6)
        else:
           sigma = coef_B
       
        if key[0] == key[3] and key[1] == key[4]:
           if not are_bonded(key[2]+1, key[5]+1, key[0], ff):
              if dist > sigma * softness:
                 energy = energy + LJ(coef_A, coef_B, dist, form)
              else:
                 if verbose:
                    print('A-self-overlap')
                 return(math.inf, math.inf)
           #else:
              #print('are bonded')
        else:
           if dist > sigma * softness:
                 energy = energy + LJ(coef_A, coef_B, dist, form)
           else:
              if verbose:
                 print('A-self-overlap')
              return(math.inf, math.inf)  
        e_pot = e_pot+ coulomb(charge_A, charge_B, dist, eps)
        
    return(energy, e_pot)
