#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                             METROPOLIS-MONTE-CARLO TOOL FOR GENERATING POLYMER CHAINS, MELTS & SOLVATED POLYMER SYSTEMS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Some random edit
import random
import math
import argparse
import itertools
import numpy as np
import scipy.optimize as opt
from numpy import sqrt, pi, cos, sin, dot, cross, arccos, degrees
from numpy import float as nfl
from numpy.linalg import norm

#======================================================================================================================================================================
#    IV                                                         DEFINITIONS & CONSTANTS
#======================================================================================================================================================================

# The force-field is globally available for efficency
global ff

# All physical constants are globally available 
global kb
kb = 1.38964852 * 10**(-23.0) # J/K

#======================================================================================================================================================================
#    V                                                       GENERAL PURPOSE FUNCTIONS
#======================================================================================================================================================================

def read_input(name):
    ff = {}
    atoms, bonds, angles, dih, constraints, virtual_sitsn = [], [], [], [], [], []
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
    ff.update({'atoms':atoms, 'bonds':bonds, 'angles':angles, 'constraints':constraints, 'dih':dih, 'nonbond_params':nonbond_params})  
    return(ff)              

def convert_constraints(STATUS):
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
    return(None)

def write_gro_file(data, name, n):
    out_file = open(name, 'w')
    out_file.write('Monte Carlo generated PS in THF'+'\n')
    out_file.write('{:>3s}{:<8d}{}'.format('',n,'\n'))
    count = 1
    atoms = ['B','R','R','R']
    res_num = 1
    index = 0
    for coord in data:
        #print(count)
        out_file.write('{:>5d}{:<5s}{:>5s}{:5d}{:8.3F}{:8.3F}{:8.3F}{}'.format(res_num, 'STYR', atoms[index], count, coord[0], coord[1], coord[2],'\n'))
        index = index + 1
        if count % 4 is 0:
            #print('--------------')
            res_num = res_num + 1
            index = 0
        count = count + 1
    out_file.write('{:>2s}{:<.5F} {:<.5F} {:<.5F}'.format('',2.000000, 2.00000, 2.00000))
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
                  coordinates[count - 3] = point
                  count = count + 1
               else:
                  count = count +1 
    return(coordinates)

#======================================================================================================================================================================
#   VI                                                       GEOMETRICAL FUNCTIONS
#======================================================================================================================================================================

def u_vect(vect):
    return(vect/norm(vect))

def angle(A, B, C):
    v1 = B - A
    v2 = B - C
    return(degrees(arccos(np.clip(dot(u_vect(v1), u_vect(v2)), -1.0, 1.0))))

def dih(A, B, C, D):
    r1 = A - B
    r2 = B - C
    r3 = C - D
    n1 = cross(r1, r2)
    n2 = cross(r2, r3)
    return(degrees(arccos(np.clip(dot(u_vect(n1), u_vect(n2)), -1.0, 1.0))))



def norm_sphere():
    v_sphere = np.random.normal(0.0, 1, (5000,3))
    return(np.array([ u_vect(vect) for vect in v_sphere]))

     
#======================================================================================================================================================================
#   VII                                                   GROMACS POTENTIAL DEFINITIONS 
#======================================================================================================================================================================


def pot_I( val ,k0, ref):
    return(0.5 * k0 * (val-ref)**2.0)

def LJ(sig, eps, r):
    return(4.0 * eps * ( (sig/r)**12.0 - (sig/r)**6.0))

def proper_dih(dih, ref, k0,n):
    return( k0 * (1 + cos(n * phi - ref)))

def legal(term, traj):
    status_A = all( [index <= len(traj) for index in term['pairs']])
    if status_A:
       coords = [traj[i - 1] for i in term['pairs']]
       status_B = all([any(coord != np.array([0,0,0])) for coord in coords])
       return(status_B)
    else:
       return(status_A)

def bonded_pot(traj):
    bond_pairs = [ norm(traj[(term['pairs'][0] - 1)] - traj[(term['pairs'][1] - 1)]) for term in ff['bonds'] if legal(term, traj)]
    return(sum([pot_I(dist, term['k0'], term['ref']) for term, dist in zip(ff['bonds'], bond_pairs)])) 

def angle_pot(traj):
    angles =  [ angle(traj[(term['pairs'][0] - 1)], traj[(term['pairs'][1] - 1)], traj[(term['pairs'][2] - 1)]) for term in ff['angles'] if legal(term, traj) ]
    return(sum([pot_I(ang, term['k0'], term['ref']) for term, ang in zip(ff['angles'], angles)])) 

def dihedral_pot(traj):
    dih_ang = [dih(traj[(t['pairs'][0] - 1)], traj[(t['pairs'][1] - 1)], traj[(t['pairs'][2] - 1)], traj[(t['pairs'][3] -1)]) for t in ff['dih'] if legal(t, traj)]
    return(sum([pot_I(ang, term['k0'], term['ref']) for term, ang in zip(ff['dih'], dih_ang)]))

# Small note here:
# Potentially we don't really wanna do dynamics with this thing. Thus we can recycle the pair list and all the other stuff since we don't move it. 
# The only thing might be that we want an in-house energy minimization after the super-CG back-transformation

def Vdw_pot(traj):
    energy=0
    for i, pos_A in enumerate(traj):
       for j, pos_B in enumerate(traj):
           dist = norm(pos_A - pos_B)
           atom_A, atom_B = ff['atoms'][i]['typ'], ff['atoms'][j]['typ']
   #        print(atom_A, atom_B  )
           epsilon = ff['nonbond_params'][(atom_A, atom_B)]['epsilon']
           sigma = ff['nonbond_params'][(atom_A, atom_B)]['sigma']
           if dist != 0:
              energy = energy + LJ(sigma, epsilon, dist)
    return(energy)

  
#======================================================================================================================================================================
#      VIII                                                 POTENTIAL ENERGY & MINIMIZATION
#======================================================================================================================================================================

def Hamiltonion(traj, display):
    traj = traj.reshape(-1,3)
    bonded = bonded_pot(traj)
    angle = angle_pot(traj)
    dihedral = dihedral_pot(traj)
    vdw = Vdw_pot(traj)
    if display:
       for term, name in zip([bonded, angle, dihedral, vdw],['bonds', 'angle', 'dihedral', 'vdw']):
           print(name, term)
    return(bonded + angle + dihedral + vdw)

def energy_min(initial, bounds):
    #print '---------bounds and initial---------------'
    #print len(bounds) 
    #print len(initial)
    opt_coord = opt.minimize(Hamiltonion, initial, method='L-BFGS-B', bounds=bounds)
    #opt_coord = opt.minimize(Hamiltonion, initial, method='CG')
  # print Hamiltonion(opt_coord.x)
    #print '------------------------'
    return(opt_coord)

#======================================================================================================================================================================
#      IX                                                      METROPOLIS-MONTE-CARLO FUNCTIONS
#======================================================================================================================================================================

def is_overlap(new_point, traj, tol, nexcl=1):
 #   print('---------- output is overlap ---------')
    n = len(traj) - nexcl
    distances = [ point - new_point for point in traj[0:n]]
  #  print(norm(distances,axis=1))
  #  print('-----------end is overlap -----------')
    return( any(norm(distances, axis=1) <  tol))

def take_step(vectors, step_length, item):
    index = random.randint(0, len(vectors) - 1)
    new_item = item + vectors[index] * step_length
    return(new_item, index)

def selv_av_random_step(traj, step_length, size, nexcl):
        #print('-> take random step ')
        subcount =0
        n=len(traj) - 1
        
        while True:
              vector_bundel = norm_sphere()
              new_coord, index = take_step(vector_bundel, step_length, traj[len(traj) - 1])
             # print('from slev_av_random')
             # print(traj)
              if n != 0:
                 if not(is_overlap(new_coord, traj, size, nexcl)):
                    break
                 elif subcount < 5000:
                    vector_bundel = np.delete(vector_bundel, index, axis=0)
                    subcount = subcount + 1
                 else:
                    print('WARNING in SELF-AV-RANDOM')
                    break
              else:
                break
        return(new_coord)

def geometrical_center(coord):
    return(sum(coord)/float(len(coord)))

def determine_step_legnth(coords, bb_indices):
       bb_coord = [coords[i] for i in bb_indices]
       step_length =  norm(sum(bb_coord))
       g0 = geometrical_center(coords)
       max_dist_from_center = max([ norm(g0 - point) for point in coords])
       size =  max_dist_from_center + 0.45
       nexcl = math.ceil(size/step_length)
       print("step:",step_length)
       print("size:", size)
       print("Will exclude",nexcl, "interactions.")
       return(step_length, size, nexcl)

#======================================================================================================================================================================
# RESOLUTION TRANSFORMATION
#======================================================================================================================================================================

def add_particels(traj, new_point, n_atoms, distances):
  #  print('-> adding particles ')
 #   print(new_point)
    #bounds = np.c_[traj.ravel(), traj.ravel()].tolist() + [[None, None]] * (3 * n_atoms)
    directions = np.random.normal(0.0, 1.0, (n_atoms,3))
    for i, direct in enumerate(directions):
           atom = new_point + u_vect(direct) * distances[i]
           traj = np.append(traj, atom)  
    #print("Froma add particles\n")
    #print(traj)
    #traj = energy_min(traj, bounds).x.reshape(-1,3)
    return(traj.reshape(-1,3))   
  
def accaptable(E, temp, prev_E):
    if E < prev_E:
       return(True)
    else:
       N = np.random.normal(0,1,(1,1))
       F = np.exp(-10**3 * (E-prev_E)/(kb*temp))
    #   F = np.exp( -(E-prev_E)/(kb*temp))
       if N < F:
          return(True)
       else:
          return(False)

def metropolis_monte_carlo(n_chains, n_repeat, conf, l_box, temp, bb_indices):
    print('\n++++++++++ Starting Monte-Carlo Program +++++++++\n')
    n_atoms = len(conf) 
    step_length, size, nexcl = determine_step_legnth(conf, bb_indices)
    traj = conf  #np.empty(0)
    cg_traj = np.array([[0.0,0.0,0.0]])
    #cg_traj = np.append(cg_traj,np.array([np.random.normal(0,1,3)]))
    count = 0                               
    prev_E = 0.0
    rejected = 0                            
    while count < n_repeat -1:                 
       #print('---', count)                  
       while True:                         
          #print(cg_traj) 
          new_cg   = selv_av_random_step(cg_traj.reshape(-1,3), step_length, size, nexcl)
          new_traj = add_particels(traj.reshape(-1,3), new_cg, n_atoms, [0])
          new_points = new_traj[-1::-(n_atoms+1)]
          new_cg  = geometrical_center(new_points)
          total_E  = Hamiltonion(new_traj, False)/len(new_traj)
          if accaptable(total_E, temp, prev_E):
           # print('accapted')
           # print(total_E * len(new_traj))
            prev_E = total_E
            traj = new_traj
            cg_traj = np.append(cg_traj, new_traj)
            count = count + 1
            #print(traj)
            #write_gro_file(traj.reshape(-1,3), str(count), len(traj))
            break
          else:
            print('rejected')
            rejected = rejected + 1
    print('++++++++++++++++ RESULTS FORM MONTE-CARLO PROGRAM ++++++++++++++\n')
    print('Total Energy:', Hamiltonion(new_traj, True))
    print('Radius of Gyration:', radius_of_gyr(traj))
    print('Number of rejected MC steps:', rejected)
    return(traj)

#======================================================================================================================================================================
#   X                                                              SOME USEFUL ANALYSIS FUNCTIONS
#======================================================================================================================================================================

def radius_of_gyr(traj):
    N = len(traj)
    diff=np.zeros((N**2))
    count=0
    for i in traj:
        for j in traj:
            diff[count]=dot((i - j),(i-j))
            count = count + 1
    Rg= 1/np.float(N)**2 * sum(diff)
    return(np.float(sqrt(Rg)))

def average_end_end_dist(traj):
    return(0)

def PDI():
    return(0)

#======================================================================================================================================================================
#   XI                                                                 MAIN & INPUT ARGUMENTS
#======================================================================================================================================================================

    
def build_system(topfile, conv, structure_file, n_chains, n_mon, box_vect, temp):
    global ff
    ff = read_input(topfile)
    convert_constraints(conv)
 #   print(ff['bonds'])
    conf = read_conf_file(structure_file, 'gro')
 #   print("Read in conf:",conf)
    traj = metropolis_monte_carlo(n_chains, n_mon, conf, box_vect, temp, [0])
    write_gro_file(traj,'out.gro',len(traj))
    return(None)

#============================================================================= END =====================================================================================
