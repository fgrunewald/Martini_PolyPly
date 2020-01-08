import math
import argparse
import numpy as np
import scipy.optimize as opt
import scipy.spatial
import scipy 
from numpy import sqrt, pi, cos, sin, dot, cross, arccos, degrees
from numpy import float as nfl
from numpy.linalg import norm
from string import digits
from polyply.structure_tool.analysis_funtions import *
from polyply.structure_tool.geometrical_functions import *
from polyply.structure_tool.force_field_tools import *
from multiprocessing import Pool
import scipy.spatial
import networkx as nx
from tqdm import tqdm

global kBa
kb = 1.38964852 * 10**(-23.0) *10**-3.0 # kJ/K


def construct_bonded_exclusions(ff):
    bonded_lists={}
    for molecule, params in ff.items():
     
      if molecule != 'nonbond_params' and molecule != 'nonbond_14pairs' and len(params['bonds']) != 0:
 #       print(molecule)
        mol_graph = construct_mol_graph(params['bonds'])
        bonded_list = {}

        for atom in params['atoms']:
            bonded_list.update({atom['n']: neighborhood(mol_graph, atom['n'], params['nexcl'])})
        bonded_lists.update({molecule:bonded_list})
      else:
        bonded_lists.update({molecule:[]})
    for molecule, bond_list in bonded_lists.items():
       params = ff[molecule]
       params.update({'bond_excl':bond_list})
       ff.update({molecule:params})
    return(ff)
  
def construct_mol_graph(bonds):
    G = nx.Graph()
    edges = [ (entry['pairs'][0], entry['pairs'][1]) for entry in bonds]
    G.add_edges_from(edges)
    return(G)

def neighborhood(G, node, n):
#     Adobted from: https://stackoverflow.com/questions/22742754/finding-the-n-degree-neighborhood-of-a-node    
    path_lengths = nx.single_source_dijkstra_path_length(G, node)
    neighbours=[ node for node, length in path_lengths.items() if length <= n]
    return(neighbours)


def read_top(name):
    itpfiles =  []
    system, ff = {}, {}
    empty = 0
    section = 'Random'
    print()
    print('{0:+^120}'.format(' reading force-field '))
    with open(name) as f:
         lines = f.readlines()
         for line in lines:
           if not line[0] in ';' and not len(line.replace('\n', '').split()) == 0:
             if len(line.replace('\n', '').split()) == 0:
                empty = empty +  1
             if is_section_head(line):
                section = line.replace('\n', '').replace(']','').replace('[','').strip()
                   #print(section)
             elif section in '[ molecules ]':
                  name, n_mol = strip_comments(line)[0:2] 
                  system.update({name: int(n_mol)})
             elif section in '[ System ]':
                  q=1
                  #print('Reading in', line)  
             elif any([ word in '#include' for word in line.split()]):
                  itpfiles += [line.replace('\n', ' ').replace('\"', ' ').split()[1]]
                  #print(line)

    for itp in itpfiles:
        parameters =  read_itp(itp)
        for key, items in parameters.items():
            if key != 'nonbond_params' and key != 'nonbond_14pairs':
               ff.update({key:items})
            else:
               ff.update({key:items[key]})
 
    ff = convert_constraints(ff)
    ff = construct_bonded_exclusions(ff)
    print('finished reading force-field')
    return(ff, system)

def strip_comments(line):
    line = line.replace('\n', ' ').split()
    clean_line = []
    count = 0
    word = 'random'
    while True:
          if count < len(line):
   #          print(count)
  #           print(len(line))
             word = line[count]
             #print(word)
             if word not in '[ ; ]' and '#' not in word:
                #print(word)
                word = line[count]
                clean_line += [word]
                count = count + 1 
             else:
               break
          else:
            break
    return(clean_line)

def is_section_head(line):

    for word in line.split():
       for char in list(word):
           if char not in '[ [ ]':
              return False
           else:
              return True

def read_itp(name):
    not_implemented={'bonds':[3,4,5,6,7,8,9,10],'angles':[3,4,5,6,7,8],'dihedrals':[8],'restraints':[],'virtual_sites':[1,2,3,4]}
    print('Reading',name)
    parameters = {}
    molecules, atoms, bonds, angles, dih, constraints, virtual_sitsn, pairs = [], [], [], [], [], [], [], []
    nonbond_params, nonbond_14pairs = {}, {}
    nexcl=0
    section = 'random'
    prev_section = 'random'
    empty=0
    line_count = 0
    with open(name) as f:
         lines=f.readlines()
         for line in lines:
           line_count = line_count + 1
           #print(line)
           if len(line.replace('\n', '').split()) == 0:
              empty = empty +  1
           elif not line[0] in ';' and "#" not in line[0]:
                try:
                     if is_section_head(line):
                        if section != 'random' and section != 'atomtypes' and section != 'defaults':
                           #print(section)
                           parameters.update({molecule_type:{'nexcl':nfl(nexcl), 'atoms':atoms, 'bonds':bonds, 'angles':angles, 'constraints':constraints, 'dih':dih,'pairs':pairs,'nonbond_params':nonbond_params,'nonbond_14pairs':nonbond_14pairs}})
                        section = line.replace('\n', '').replace(']','').replace('[','').strip()
                        print(section)
                     
                     elif section in '[ moleculetype ]':
                          name, nexcl = strip_comments(line)[0:2]
                          molecules.append({'name':name,'nexcl':nfl(nexcl)})
                          molecule_type = name
                
                     elif section in '[ atomtypes ]':
                          atom = line.replace('\n', '').split()[0]
                          sigma, epsilon = strip_comments(line)[-2:]
                          nonbond_params.update({(atom, atom): {'sigma':nfl(sigma), 'epsilon':nfl(epsilon)}})
                           
                     elif section in '[ atoms ]':
                         n, typ, resnr, res, atom, cgnr, charge = strip_comments(line)[:7]
                         atoms.append({'n': int(n), 'typ':typ ,'atom':atom, 'charge':nfl(charge), 'resname':res})
                     
                     elif section in '[ nonbond_params ]':
                         atom1, atom2, f, sigma, epsilon = strip_comments(line)
                         nonbond_params.update({(atom1, atom2): {'sigma':nfl(sigma), 'epsilon':nfl(epsilon)}})
                         molecule_type = 'nonbond_params'
                    
                     elif section in '[ virtual_sitesn ]':
                         clean_line = strip_comments(line)
                         
                         # this is a workaround for ignoring VS, we still can't construct them
                         atom1 = clean_line[0]
                         for atom2 in clean_line[1:]:
                             bonds.append({'pairs':[int(atom1),int(atom2)], 'k0':nfl(10000), 'ref':nfl(0.43), 'f':nfl(1)})
                
                     elif section in '[ pairtypes ]':
                         atom1, atom2, f, sigma, epsilon = strip_comments(line)
                         nonbond_14pairs.update({(atom1, atom2): {'sigma':nfl(sigma), 'epsilon':nfl(epsilon)}})
                         molecule_type = 'nonbond_14pairs'                        
 
                     elif section in '[ pairs ]':
                         A, B, f = strip_comments(line)
                         pairs.append({'pairs':(int(A),int(B)), 'f':nfl(f)})
                
                     elif section in '[ bonds ]':
                         A, B, f, ref, k0 = strip_comments(line)
                         bonds.append({'pairs':[int(A),int(B)], 'k0':nfl(k0), 'ref':nfl(ref), 'f':nfl(f)})
                     
                     elif section in '[ angles ]':
                         A, B, C, f, ref, k0 = strip_comments(line)
                         angles.append({'pairs': [int(A), int(B), int(C)], 'k0':nfl(k0), 'ref':nfl(ref),'f':nfl(f)})
                     
                     elif section in '[ dihedrals ]':
                         if line.split()[4] == str(1) or line.split()[4] == str(9):
                            A, B, C, D, f, ref, k0, n = strip_comments(line)
                            dih.append({'pairs':[int(A),int(B),int(C), int(D)], 'k0':nfl(k0), 'f':nfl(f), 'n':nfl(n), 'ref':nfl(ref)})
                         elif line.split()[4] == str(2) or line.split()[4] == str(10):
                            A, B, C, D, f, ref, k0 = strip_comments(line)
                            dih.append({'pairs':[int(A),int(B),int(C), int(D)], 'k0':nfl(k0), 'f':nfl(f), 'n':0, 'ref':nfl(ref)})
                         elif line.split()[4] == str(3) or line.split()[4] == str(5) or line.split[4] == str(11):
                            A, B, C, D, f = strip_comments(line)[0:5]
                            Clist = [ nfl(value) for value in strip_comments(line)[5:]]
                            dih.append({'pairs':[int(A),int(B),int(C), int(D)], 'clist':Clist,'f':nfl(f)})
                         elif line.split()[4] == str(4):
                            A, B, C, D, f, ref, k0, n = strip_comments(line)
                            dih.append({'pairs':[int(A),int(B),int(C), int(D)], 'ref':nfl(ref),'k0':nfl(k0),'n':nfl(n),'f':nfl(f)})
                                      
                     elif section in '[ constraints ]':
                         A, B, f, ref = strip_comments(line)
                         constraints.append({'pairs':[A, B], 'f':nfl(f), 'ref':nfl(ref)})
                     elif section in '[ defaults ]':
                         LJ = int(line.replace('\n', '').split()[1])
                         if LJ == 1:
                            form='C6C12'
                         else:
                            form='sigeps'
                         nonbond_params.update({'functype': form})
                except ValueError:
                       print('{0:+^80}'.format(' Fatal Error '))
                       print('Too few or too many paramters found in file ',name,'on line ',line_count,'.')
                       print('The relevant topolgy section is: ',section)
                       print('The follwing function types of this section are currently not implemented:')
                       print(not_implemented[section])
                       print('Please make sure to provide the same number and order of paramteres as outlined \n')
                       print('in the GROMACS Manual 2016.3 page 136-140.')
                       print('If you are sure your format is correct please submit a bug report on github.')
                       exit()
    if section != 'random' and section != 'atomtypes' and section != 'defaults':
       parameters.update({molecule_type:{'nexcl':nfl(nexcl), 'atoms':atoms, 'bonds':bonds, 'angles':angles, 'constraints':constraints, 'dih':dih,'pairs':pairs,'nonbond_params':nonbond_params,'nonbond_14pairs':nonbond_14pairs}})
    return(parameters)

def convert_constraints(ff):
    for molecule in ff:
     if molecule != 'nonbond_params' and molecule != 'nonbond_14pairs':
      if len(ff[molecule]['constraints']) != 0:
        print('converting constraints to bonds')
        new_bonds = []
        #print(len(ff[molecule]['bonds']))
        new_bonds = [ {'pairs':[int(term['pairs'][0]), int(term['pairs'][1])], 'ref':term['ref'], 'k0': nfl(9000.0),'f':1} for term in ff[molecule]['constraints'] ]
        new_bonds =  ff[molecule]['bonds'] + new_bonds
        ff[molecule].update({'bonds':new_bonds})
    return(ff)

def write_gro_file(data, name, ff, box):
    n = sum([ len(coords) for resname, list_of_coords in data.items() for coords in list_of_coords ])
    out_file = open(name, 'w')
    out_file.write('Monte Carlo generated PEO'+'\n')
    out_file.write('{:>3.3s}{:<8d}{}'.format('',n,'\n'))
    count = 0
    resnum = 1
  
    for mol_name, mols in data.items():
       for coords in mols:
         resnum = resnum + 1
         for index, line in enumerate(coords):
             atomtype = ff[mol_name]['atoms'][index]['atom']
             resname  = ff[mol_name]['atoms'][index]['resname']
             count = count + 1
             out_file.write('{:>5d}{:<5.5s}{:>5.5s}{:5d}{:8.3F}{:8.3F}{:8.3F}{}'.format(resnum, resname, atomtype, count, line[0], line[1], line[2],'\n'))
    out_file.write('{:>2s}{:<.5F} {:<.5F} {:<.5F}'.format('',float(box[0]), float(box[1]), float(box[2])))
    out_file.close()
    return(None)

def pot_I( val ,k0, ref):
    return(0.5 * k0 * (val-ref)**2.0)   

def proper_dih(dih, k0, ref, n, f):
    if f == 1 or f == 9:
       return( k0 * (1 + cos(n * np.radians(dih) - np.radians(ref))))
    elif f == 2:
       return(0.5 * k0 * (dih - ref)**2.0 )

def bond_pot(val,k0,ref,f):

   
  
    if f == 1:
       val = 0.5 * k0 * (val-ref)**2.0
    elif f == 2:
       val = 0.25 * k0 * (val**2.0-ref**2.0)**2.0

    #if val == None:
    #   print(val,k0,ref,f) 

    return(val)

def ang_pot(val,k0,ref,f):
    if f == 1:
       return(0.5 * k0 * (val-ref)**2.0)
    elif f == 2:
       #print('go here') 
       return(0.5 * k0 * (cos(val)-cos(ref))**2.0)
    elif f == 10:
       return(0.5 * k0 * (cos(val)-cos(ref))**2.0/sin(val)**2.0)

def legal(term, traj, restart):
    status_A = all( [index <= len(traj) for index in term['pairs']])
    status_B = all( [index >= restart for index in term['pairs']])
    
    if status_A and status_B:
       return True
    else:
       return False 

def bonded_pot(ff, traj, restart):
    #print(restart)
    bond_pairs = [ norm(traj[(term['pairs'][0] - 1)] - traj[(term['pairs'][1] - 1)]) for term in ff['bonds'] if legal(term, traj, restart)]   
    #print([(dist, term['k0'], term['ref'],term['f']) for term, dist in zip(ff['bonds'], bond_pairs)])

    return(sum([bond_pot(dist, term['k0'], term['ref'],term['f']) for term, dist in zip(ff['bonds'], bond_pairs)]))

def angle_pot(ff, traj, restart):
    angles =  [ angle(traj[(term['pairs'][0] - 1)], traj[(term['pairs'][1] - 1)], traj[(term['pairs'][2] - 1)]) for term in ff['angles'] if legal(term, traj, restart) ]
    return(sum([ang_pot(np.radians(ang),term['k0'], np.radians(term['ref']),term['f']) for term, ang in zip(ff['angles'], angles)]))

def dihedral_pot(ff, traj, restart):
    dih_ang = [dih(traj[(t['pairs'][0] - 1)], traj[(t['pairs'][1] - 1)], traj[(t['pairs'][2] - 1)], traj[(t['pairs'][3] -1)]) for t in ff['dih'] if legal(t, traj, restart)]
    return(sum([proper_dih(ang, term['k0'], term['ref'], term['n'], term['f']) for term, ang in zip(ff['dih'], dih_ang)]))


def are_bonded(atom_A, atom_B, molecule, ff):
    return(any([ atom_B == atom for atom in ff[molecule]['bond_excl'][atom_A]]))

def are_14(atom_A, atom_B, molecule, ff):
    #print(ff[molecule]['nonbond_14pairs'])
    return(any([ set((atom_A, atom_B)) == set(ref['pairs']) for ref in ff[molecule]['pairs']]))

def coulomb(ca, cb, dist,eps):
    return(138.935458/(eps)*(ca*cb)/dist)

def LJ(A, B, r, form):
    if form == 'C6C12':
       return( A/r**12.0 - B/r**6.0 )
    elif form == 'sigeps':
       #print('go here', r, A,B) 
       return(4.0 * A * ( (B/r)**12.0 - (B/r)**6.0) )

def lookup_interaction_parameters(ff, atom_A, atom_B, LJ_form, key):

      try:
          coef_A = ff[key][(atom_A, atom_B)]['epsilon']
          coef_B = ff[key][(atom_A, atom_B)]['sigma']
      except KeyError:
          coef_A = ff[key][(atom_B, atom_A)]['epsilon']
          coef_B = ff[key][(atom_B, atom_A)]['sigma']

      if LJ_form == 'C6C12':
         if coef_A != 0:
            sigma = (coef_A/coef_B)**(1/6)
         else:
            sigma = 0
      else:
         sigma = coef_B

      return(coef_A, coef_B, sigma)

def nonbonded_potential(dist_mat, ff, softness, eps, form, verbose, offset):
    LJ_energy = 0
    COUL_energy=0
   # verbose=True
   # print(len(dist_mat))
   
    for key, dist in dist_mat.items():
        #print('go')
        atom_A, atom_B = ff[key[0]]['atoms'][key[2]]['typ'], ff[key[3]]['atoms'][key[5]]['typ']
        
        if key[2] < offset and key[5] < offset:
           #print('go here') 
           continue 

        charge_A, charge_B = ff[key[0]]['atoms'][key[2]]['charge'], ff[key[3]]['atoms'][key[5]]['charge']
        
        if key[0] == key[3] and key[1] == key[4]:
         #  print('go')
           if not are_bonded(key[2]+1, key[5]+1, key[0], ff) and not are_14(key[2]+1, key[5]+1, key[0], ff):
              coef_A, coef_B, sigma = lookup_interaction_parameters(ff, atom_A, atom_B, form, 'nonbond_params')

              if dist > sigma * softness:
                 LJ_energy = LJ_energy + LJ(coef_A, coef_B, dist, form)
               #  print( LJ(coef_A, coef_B, dist, form))
                 COUL_energy = COUL_energy + coulomb(charge_A, charge_B, dist, eps)
              else:
                 if verbose:
                    print('A-self-overlap')
                 return(math.inf, math.inf)

           elif are_14(key[2]+1, key[5]+1, key[0], ff):
                coef_A, coef_B, sigma = lookup_interaction_parameters(ff, atom_A, atom_B, form, 'nonbond_14pairs')
              #  print(dist, sigma*softness)
                if dist > sigma * softness: 
                   LJ_energy = LJ_energy + LJ(coef_A, coef_B, dist, form) 
                   COUL_energy = COUL_energy + coulomb(charge_A, charge_B, dist, eps)
                else:
                   return(math.inf, math.inf)
              

        else:
           coef_A, coef_B, sigma = lookup_interaction_parameters(ff, atom_A, atom_B, form, 'nonbond_params')
           if dist > sigma * softness:
               LJ_energy = LJ_energy + LJ(coef_A, coef_B, dist, form)
               print(LJ_energy)
               COUL_energy = COUL_energy + coulomb(charge_A, charge_B, dist, eps)

           else:
              if verbose:
                 print('A-self-overlap')
              return(math.inf, math.inf)  
        
    return(LJ_energy, COUL_energy)
