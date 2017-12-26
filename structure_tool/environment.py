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


def create_environment(options):
    env_type, sol, box, spacing, lipid_type = options 

    if env_type in '[ bulk ]':
       supp_traj = make_box_w_solvent(options)

    elif env_type in '[ bilayer ]':
        area_per_lipid = 0.8 # this should be made variable but on average 0.8 should be fine
        z_dim = box[2]
        xy_dim = spacing / area_per_lipid
        n_lipids = int(xy_dim**2.0)
        supp_traj = gen_bilayer(lipid_type, n_lipids, xy_dim, xy_dim, z_dim)
      
    return(supp_traj)

def gen_bilayer(lipid_type, n_lipids, x_dim, y_dim, z_dim):

    '''
    This module returns a ready to go bilayer made by insane. 
    We use suprocess because insane is currently in python 2. 
    '''
    variable_lipid = lipid_type + ":" + str(n_lipids - 1)
    
    insane_command = "insane -o bilayer.gro -p topol.top -x " + str(x_dim) + " -y " + str(y_dim) + " -z " + str(z_dim) + " -l " + str(variable_lipid) + " -l DOPE:1 -salt 0 -sol W"
    print(insane_command)
    process = subprocess.Popen(insane_command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    bilayer_coords = read_conf_file('bilayer.gro', ".gro")
    return(bilayer_coords)

def solvate(options):
    solvated_box=0
    return(solvated_box)
