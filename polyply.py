import argparse
from Martini_PolyPly.itp_tool.itp_I import *
from Martini_PolyPly.structure_tool.mc_poly_growth import *

parser = argparse.ArgumentParser(description=' ')
parser.add_argument('-itp'   , metavar = 'itps of blocks'  , dest = 'itpfiles', type = str   , help = 'name of the monomer .itp files', nargs = '*')
parser.add_argument('-n_mon' , metavar = 'length of blocks', dest = 'mon'     , type = int   , help = 'number of monomers per blocks', nargs = '*')
parser.add_argument('-name'  , metavar = 'name of molecule', dest = 'name'    , type = str   , help = 'name of the new polymer molecule', default='polymer')
parser.add_argument('-o'     , metavar = 'name of outfile' , dest = 'outfile' , type = str   , help = 'name of the new .itp file.' )
parser.add_argument('-p'     , metavar = 'topology file'   , dest = 'topfile' , type = str   , help = 'name of system topology file')
parser.add_argument('-c'     , metavar = 'structure file'  , dest = 'grofile' , type = str   , help = 'name of monomer structure file')
parser.add_argument('-T'     , metavar = 'temperature'     , dest = 'temp'    , type = float , help = 'temperature of system', default=298.0) 
parser.add_argument('-sys'   , metavar = 'system'          , dest = 'system'  , type = str   , help = 'type of system to create' , default='vac')
parser.add_argument('-v'     , metavar = 'verbose'         , dest = 'v'       , type = bool  , help = 'be loud and noisy', default=False)
parser.add_argument('-conv'  , metavar = 'convert-constraints', dest = 'conv' , type = bool  , help = 'convert constraints to bonds for minimization', default=True)
#parser.add_argument('-excl'  , metavar = 'n_excluded'      , dest = 'nexcl'   , type = int   , help = 'number of excluded interactions', default=3)
args = parser.parse_args()

def main():
   if not args.itpfiles == None:
      itp_tool(args.itpfiles, args.mon, args.outfile, args.name)
   elif not args.system == None:
      build_system(args.topfile, args.conv, args.grofile, 1, args.mon[0], np.array([5.0,5.0,5.0]), args.temp)
   else:
      print('Please specify either -itp or -sys option.')
   return(None)

main()
