import random
import argparse
import itertools
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import Martini_PolyPly.itp_tool.itp_I
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
from numpy import sqrt, pi, cos, sin, dot, cross, arccos, degrees
from numpy import float as nfl
from numpy.linalg import norm


parser = argparse.ArgumentParser(description=' ')
parser.add_argument('-itp'   , metavar = 'itps of blocks'  , dest = 'itpfiles', type = str   , help = 'name of the monomer .itp files', nargs = '*')
parser.add_argument('-n_mon' , metavar = 'length of blocks', dest = 'mon'     , type = int   , help = 'number of monomers per blocks', nargs = '*')
parser.add_argument('-o'     , metavar = 'name of outfile' , dest = 'outfile' , type = str   , help = 'name of the new .itp file.' )
parser.add_argument('-p'     , metavar = 'topology file'   , dest = 'topfile' , type = str   , help = 'name of system topology file')
parser.add_argument('-c'     , metavar = 'structure file'  , dest = 'grofile' , type = str   , help = 'name of monomer structure file')
parser.add_argument('-T'     , metavar = 'temperature'     , dest = 'temp'    , type = float , help = 'temperature of system', default=298.0) 
parser.add_argument('-sys'   , metavar = 'system'          , dest = 'system'  , type = str   , help = 'type of system to create' , default='vac')
parser.add_argument('-v'     , metavar = 'verbose'         , dest = 'v'       , type = bool  , help = 'be loud and noisy', default=False)
parser.add_argument('-excl'  , metavar = 'n_excluded'      , dest = 'nexcl'   , type = int   , help = 'number of excluded interactions', default=3)
args = parser.parse_args()

def main():
   if not args.itpfiles == None:
      itp_tool()
   elif not args.sys == None:
      build_system()
   else:
      print('Please specify either -itp or -sys option.')
   return(None)

main()
