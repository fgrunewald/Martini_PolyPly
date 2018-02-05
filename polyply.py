import argparse
import os
from Martini_PolyPly.itp_tool.itp_I import *
#from Martini_PolyPly.monomer_itps.resolve_name import *
from Martini_PolyPly.structure_tool.mc_poly_growth import *
from Martini_PolyPly.structure_tool.analysis_funtions import *
from Martini_PolyPly.structure_tool.geometrical_functions import *
from Martini_PolyPly.structure_tool.force_field_tools import *


parser = argparse.ArgumentParser(description=' ')
parser.add_argument('-itp'     , metavar = 'itps of blocks'     , dest = 'itpfiles', type = str   , help = 'name of the monomer .itp files', nargs = '*')
parser.add_argument('-n_mon'   , metavar = 'length of blocks'   , dest = 'mon'     , type = int   , help = 'number of monomers per blocks', nargs = '*')
parser.add_argument('-polymer' , metavar = 'itp from monomers'  , dest = 'r_itp_name', type = str   , help = 'use itps in default repository', nargs='*')
parser.add_argument('-name'    , metavar = 'name of molecule'   , dest = 'name'    , type = str   , help = 'name of the new polymer molecule', default='polymer')
parser.add_argument('-o'       , metavar = 'name of outfile'    , dest = 'outfile' , type = str   , help = 'name of the output file', default='out.gro')
parser.add_argument('-p'       , metavar = 'topology file'      , dest = 'topfile' , type = str   , help = 'name of system topology file')
parser.add_argument('-T'       , metavar = 'temperature'        , dest = 'temp'    , type = float , help = 'temperature of system', default=298.0) 
parser.add_argument('-v'       , metavar = 'verbose'            , dest = 'v'       , type = bool  , help = 'be loud and noisy')
parser.add_argument('-env'     , metavar = 'environmnet'        , dest = 'env'     , type = str   , help = 'type of environment to add to chain', default=None)
parser.add_argument('-sys'     , metavar = 'system'             , dest = 'sys'     , type = str   , help = 'trajectory of environment',default=None)   
parser.add_argument('-maxsteps', metavar = 'max MC steps'       , dest = 'maxsteps', type = int   , help = 'maximum number of MC steps befaure exit by error', default=5000)
parser.add_argument('-links'   , metavar = 'linkfile'           , dest = 'linkfile', type = str   , help = 'file where the bonded links are specified.', default=None)
parser.add_argument('-endgroup', metavar = 'endgroup itps'     , dest = 'endgroup', type = str   , help = 'itp files with endgroups', nargs='*' ,default=None)
parser.add_argument('-lipid'   , metavar = 'lipid type'        , dest = 'lipid'   , type = str   , help = 'lipid type of the bilayer when env is bilayer')
#parser.add_argument('-box'     , metavar = 'dimensions of box' , dest = 'box'     , type = float , help = 'boxdimensions in nm for x, y, z', nargs=3)
parser.add_argument('-spacing' , metavar = 'spacing of PEL'    , dest = 'spacing' , type = float , help = 'spaceing between PEL in a bilayer', default=0 )
parser.add_argument('-sol'     , metavar = 'solvent'           , dest = 'solvent' , type = str   , help = 'name of solvent for the polymer system', default=None)   
parser.add_argument('-lib'     , action='store_true'  , help = 'show itp library')   
args = parser.parse_args()


def resolve_name(names, path,ff):
    paths = []
    for name in names:
        full_name = path + 'monomer_itps/' + name + '.' + ff + '.itp'
        paths += [ full_name ]
    return(paths)

def show_files(path):
    files = [ f for f in os.listdir(path) if len(f.split('.')) == 3] # if os.path.isfile(f)]
    print('\nPolymer       ','ForceField')
    print('----------------------------')
    for f in files:
        print('{:<15s}{:<15s}'.format(f.split('.')[0], f.split('.')[1]))
  
def main():
 
   if not args.itpfiles == None:
      itp_tool(args.itpfiles, args.linkfile ,args.mon, args.outfile, args.name, args.endgroup)

   elif not args.r_itp_name == None:
      path = os.path.abspath(__file__).replace('polyply.py', '')
      itp_files = resolve_name(args.r_itp_name, path, 'martini')
      itp_tool(itp_files, args.linkfile ,args.mon, args.outfile, args.name, args.endgroup)

   elif not args.env == None:
      env_options = (args.env, args.solvent, args.lipid, args.sys)
      top_options = (args.topfile)
      mc_options =  (args.temp, args.maxsteps, args.v, args.name)
      build_system(top_options, env_options, mc_options, args.outfile)

   elif args.lib:
      path = os.path.abspath(__file__).replace('polyply.py', '')
      path = path + 'monomer_itps'           
      show_files(path)
   else:
      print('Please specify either -itp or -sys option.')

   return(None)

main()
