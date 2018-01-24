import subprocess
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

def find_central_starting_point(coordinates, sol):
    x_coords = [ atom[0] for mol in coordinates[sol] for atom in mol ]
    y_coords = [ atom[1] for mol in coordinates[sol] for atom in mol ]
    z_coords = [ atom[2] for mol in coordinates[sol] for atom in mol ]
    point_A = np.array([max(x_coords), max(y_coords), max(z_coords)])
    point_B = np.array([min(x_coords), min(y_coords), min(z_coords)])    
    diagonal = point_A - point_B
    starting_point = point_B + 0.5 * diagonal
    return(starting_point)

def import_environment(options):
    env_type, sol, lipid_type, sysfile = options 
    environment_coords, box = read_conf_file(sysfile, ".gro")
 
    if env_type in '[ sol ]':
         start = find_central_starting_point(environment_coords, sol)
         constraints = [{'type':None}]

    elif env_type in '[ bilayer ]':
         #DOPE_new = reorder_lipid(bilayer_coords['DOPE'][0])
         #environment_coords['DOPE'][0] = DOPE_new
         start = bilayer_coords[lipidtype][0][-1]
         constraints = [{'type':'dist-x-axis', 'tol':float(box[0])/2, 'ref':start[0]}, {'type':'dist-y-axis', 'tol':float(box[1])/2, 'ref':start[1]},  
                        {'type':'dist-z-axis', 'tol':0, 'ref':start[2]}]

    return(environment_coords, constraints, start, box)


def reorder_lipid(lipid_coords):
    good_lipid = [0,0,0,0,0,0,0,0,0,0,0]
    good_lipid[0:7] = lipid_coords[4:12]
    good_lipid[8] = lipid_coords[3]
    good_lipid[9] = lipid_coords[2]
    good_lipid[10] = lipid_coords[1]
    good_lipid[11] = lipid_coords[0]
    return(good_lipid)

def gen_bilayer(lipid_type, n_lipids, dim):

    '''
    This module returns a ready to go bilayer made by insane. 
    We use suprocess because insane is currently in python 2. 
    '''

    variable_lipid = lipid_type + ":" + str(n_lipids - 1)
    x_dim, y_dim, z_dim = dim
    insane_command = "insane -o bilayer.gro -p topol.top -x " + str(x_dim) + " -y " + str(y_dim) + " -z " + str(z_dim) + " -l " + str(variable_lipid) + " -l DOPE:1 -salt 0 -sol W"
    
    print(insane_command)
    
    process = subprocess.Popen(insane_command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    
    grompp_command = "gmx_mpi grompp -f min.mdp -c bilayer.gro -p topol.top"
    process = subprocess.Popen(grompp_command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    
    mdrun_command = "gmx_mpi mdrun -v"
    process = subprocess.Popen(mdrun_command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

    mdrun_command = "gmx_mpi trjconv -f confout.gro -o out.gro -pbc whole"
    process = subprocess.Popen(mdrun_command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

    bilayer_coords, box = read_conf_file('out.gro', ".gro")  
    DOPE_new = reorder_lipid(bilayer_coords['DOPE'][0])
    #print('New coordinates')
    bilayer_coords['DOPE'][0] = DOPE_new
    head = bilayer_coords['DOPE'][0][11]
    #for atom in DOPE_new:
    #    print(10*atom)

    constraints = [{'type':'dist-x-axis', 'tol':float(box[0])/2, 'ref':head[0]}, {'type':'dist-y-axis', 'tol':float(box[1])/2, 'ref':head[1]},  
                   {'type':'dist-z-axis', 'tol':0, 'ref':head[2]}]
    return(bilayer_coords, constraints, head, box)
