import argparse
import os
from Martini_PolyPly.itp_tool.itp_I import *
from Martini_PolyPly.monomer_itps.resolve_name import *
from Martini_PolyPly.structure_tool.mc_poly_growth import *
from Martini_PolyPly.structure_tool.analysis_funtions import *
from Martini_PolyPly.structure_tool.geometrical_functions import *
from Martini_PolyPly.structure_tool.force_field_tools import *


parser = argparse.ArgumentParser(description=' ')
parser.add_argument('-itp'    , metavar = 'itps of blocks'     , dest = 'itpfiles', type = str   , help = 'name of the monomer .itp files', nargs = '*')
parser.add_argument('-n_mon'  , metavar = 'length of blocks'   , dest = 'mon'     , type = int   , help = 'number of monomers per blocks', nargs = '*')
parser.add_argument('-polymer', metavar = 'itp from monomers'  , dest = 'r_itp_name', type = str   , help = 'use itps in default repository', nargs='*')
parser.add_argument('-name'   , metavar = 'name of molecule'   , dest = 'name'    , type = str   , help = 'name of the new polymer molecule', default='polymer')
parser.add_argument('-o'      , metavar = 'name of outfile'    , dest = 'outfile' , type = str   , help = 'name of the new .itp file.' )
parser.add_argument('-p'      , metavar = 'topology file'      , dest = 'topfile' , type = str   , help = 'name of system topology file')
parser.add_argument('-c'      , metavar = 'structure file'     , dest = 'grofile' , type = str   , help = 'name of monomer structure file')
parser.add_argument('-T'      , metavar = 'temperature'        , dest = 'temp'    , type = float , help = 'temperature of system', default=298.0) 
parser.add_argument('-chain'  , metavar = 'chain topology'     , dest = 'chain'   , type = str   , help = 'topology of chain' , default='linear-random')
parser.add_argument('-v'      , metavar = 'verbose'            , dest = 'v'       , type = bool  , help = 'be loud and noisy')
parser.add_argument('-conv'   , metavar = 'convert-constraints', dest = 'conv'    , type = bool  , help = 'convert constraints to bonds for minimization', default=True)
parser.add_argument('-env'    , metavar = 'environmnet'        , dest = 'env'     , type = str   , help = 'type of environment to add to chain', default='vac')
parser.add_argument('-maxsteps', metavar = 'max MC steps'       , dest = 'maxsteps', type = int   , help = 'maximum number of MC steps befaure exit by error', default=5000)
parser.add_argument('-links'  , metavar = 'linkfile'           , dest = 'linkfile', type = str   , help = 'file where the bonded links are specified.', default=None)
parser.add_argument('-endgroup', metavar = 'endgroup itps'     , dest = 'endgroup', type = str   , help = 'itp files with endgroups', nargs='*' ,default=None)
parser.add_argument('-lipid'   , metavar = 'lipid type'        , dest = 'lipid'   , type = str   , help = 'lipid type of the bilayer when env is bilayer')
parser.add_argument('-box'     , metavar = 'dimensions of box' , dest = 'box'     , type = float , help = 'boxdimensions in nm for x, y, z', nargs=3)
parser.add_argument('-spacing' , metavar = 'spacing of PEL'    , dest = 'spacing' , type = float , help = 'spaceing between PEL in a bilayer', default=0 )
parser.add_argument('-sol'     , metavar = 'solvent'           , dest = 'solvent' , type = str   , help = 'name of solvent for the polymer system', default=None)   
args = parser.parse_args()

def main():
 
   if not args.itpfiles == None:
      itp_tool(args.itpfiles, args.linkfile ,args.mon, args.outfile, args.name, args.endgroup)

   elif not args.r_itp_name == None:
      path = os.path.abspath(__file__).replace('polyply.py', '')
      itp_files = resolve_name(args.r_itp_name, path)
      itp_tool(itp_files, args.linkfile ,args.mon, args.outfile, args.name, args.endgroup)

   elif not args.chain == None:
      env_options = (args.env, args.solvent, args.box, args.spacing, args.lipid)
      top_options = (args.topfile, args.grofile, args.chain, args.conv)
      mc_options =  (args.mon[0], args.temp, args.maxsteps, args.v)
      build_system(top_options, env_options, mc_options)

   else:
      print('Please specify either -itp or -sys option.')

   return(None)

main()
