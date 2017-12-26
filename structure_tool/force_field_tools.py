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


global kB
kb = 1.38964852 * 10**(-23.0) *10**-3.0 # kJ/K

def read_itp(name):
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
              #print line.replace('\n', '')
              if any([ word in '[ [ ]' for word in line.split()]):
                 section = line.replace('\n', '').split()[1]
                 #print section
              elif section in '[ moleculetype ]':
                   name, nexcl = line.replace('\n', '').split()
                   molecules.append({'name':name,'nexcl':nfl(nexcl)})
              elif section in '[ atoms ]':
                  # print line
                  n, typ, resnr, res, atom, cgnr, charge, mass = line.replace('\n', '').split()
                  atoms.append({'n': int(n), 'typ':typ ,'atom':atom, 'charge':nfl(charge)})
              elif section in '[ nonbond_params ]':
                  atom1, atom2, sigma, epsilon = line.replace('\n', '').split()
                  nonbond_params.update({(atom1, atom2): {'sigma':nfl(sigma), 'epsilon':nfl(epsilon)}})
              elif section in '[ bonds ]':
                  # print line
                  A, B, f, ref, k0 = line.replace('\n', '').split()
                  bonds.append({'pairs':[int(A),int(B)], 'k0':nfl(k0), 'ref':nfl(ref)})
              elif section in '[ angles ]':
                  #print line
                  A, B, C, f, ref, k0 = line.replace('\n', '').split()
                  angles.append({'pairs': [int(A), int(B), int(C)], 'k0':nfl(k0), 'ref':nfl(ref)})
                  #angles.append({'pairs': [int(A), int(B), int(C)], 'k0':1000.0, 'ref':nfl(ref)})
              elif section in '[ dihedrals ]':
                  # print line
                  A, B, C, D, f, ref, k0, n = line.replace('\n', '').split()
                  dih.append({'pairs':[int(A),int(B),int(C), int(D)], 'k0':nfl(k0), 'f':nfl(f), 'n':nfl(n), 'ref':nfl(ref)})
              elif section in '[ constraints ]':
                  #print line.replace('\n', '')
                  A, B, f, ref = line.replace('\n', '').split()
                  constraints.append({'pairs':[A, B], 'f':nfl(f), 'ref':nfl(ref)})
    ff.update({'atoms':atoms, 'bonds':bonds, 'angles':angles, 'constraints':constraints, 'dih':dih, 'nonbond_params':nonbond_params, 'molecules':molecules})
    return(ff)

def convert_constraints(ff, STATUS):
    # This is not really fast but mehh it works
    # For very many molecules we should make this more efficent
    if STATUS:
       new_bonds=[]
       new_bonds = [ {'pairs':[int(term['pairs'][0]), int(term['pairs'][1])], 'ref':term['ref'], 'k0': nfl(8000.0)} for term in ff['constraints'] ]
       new_bonds =  ff['bonds'] + new_bonds
       ff.update({'bonds':new_bonds})
    else:
       print("!!!!!!!!!!! WARNING !!!!!!!!!")
       print("Constraints don't minimize that well!")
       print("Your structure might be quite distorted!")
       print("If you want to keep them constraints, use a different program.")
       exit()
    return(ff)

def write_gro_file(data, name, n):
    out_file = open(name, 'w')
    out_file.write('Monte Carlo generated PEO'+'\n')
    out_file.write('{:>3s}{:<8d}{}'.format('',n,'\n'))
    count = 1
    res_num = 1
    for coord in data:
        out_file.write('{:>5d}{:<5s}{:>5s}{:5d}{:8.3F}{:8.3F}{:8.3F}{}'.format(1, 'PEO', 'PEO', count, coord[0], coord[1], coord[2],'\n'))
        count = count + 1
    out_file.write('{:>2s}{:<.5F} {:<.5F} {:<.5F}'.format('',10.000000, 10.00000, 10.00000))
    out_file.close()
    return(None)

def write_log_file(infos, fomat_info):
    log_file = open('logfile.log', 'w')
    for line in infos:
        log_file.write(format_info.format(info))
    log_file.close()
    return(None)

def read_conf_file(filename, file_type):
    with open(filename) as f:
        lines=f.readlines()
        coordinates=np.zeros(((len(lines)-3),3))
        res_num_names = []
        count=0
        if file_type in '[ .gro, gro]':
           for line in lines:
               if count == 0:
                  title = line.replace('\n', '')
                  print(title)
                  count = count + 1
               elif 1 < count < len(lines) - 1:
                  print(line.replace('\n', '').split())
                  #res_num, res_name, atom, a_index, x, y, z, v1, v2, v3 = line.replace('\n', '').split()
                  # In principle one can also put velocities in the gro file so we should account for that at some point
                  res_num_name, atom, a_index, x, y, z = line.replace('\n', '').split()
                  point = np.array([x,y,z])
                  res_num_names += res_num_name
                  coordinates[count - 3] = point
                  count = count + 1
               else:
                  count = count +1
    return(coordinates)


def pot_I( val ,k0, ref):
    return(0.5 * k0 * (val-ref)**2.0)

def LJ(sig, eps, r):
  #  print(sig)
 #   print(eps)
    return(4.0 * eps * ( (sig/r)**12.0 - (sig/r)**6.0))
   
def proper_dih(dih, k0, ref, n):
    return( k0 * (1 + cos(n * np.radians(dih) - np.radians(ref))))

def legal(term, traj):
    status_A = all( [index <= len(traj) for index in term['pairs']])
    if status_A:
       coords = [traj[i - 1] for i in term['pairs']]
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
    #print('Dihedral')
    #print(dih_ang)
    return(sum([proper_dih(ang, term['k0'], term['ref'], term['n']) for term, ang in zip(ff['dih'], dih_ang)]))

def Vdw_pot(ff, traj):
    energy=0
    for i, pos_A in enumerate(traj):
       for j, pos_B in enumerate(traj):
           dist = norm(pos_A - pos_B)
     #      print(dist)
           atom_A, atom_B = ff['atoms'][i]['typ'], ff['atoms'][j]['typ']
   #        print(atom_A, atom_B  )
           epsilon = ff['nonbond_params'][(atom_A, atom_B)]['epsilon']
           sigma = ff['nonbond_params'][(atom_A, atom_B)]['sigma']
           if i-j > ff['molecules'][0]['nexcl']:
              energy = energy + LJ(sigma, epsilon, dist)
        #      print('Positions')
        #      print(pos_A)
        #      print(pos_B)
        #      print(dist)
        #      print(LJ(sigma, epsilon, dist))
        #      print('-----------\n')
             #  print(i,' ' , j)

    return(energy)
