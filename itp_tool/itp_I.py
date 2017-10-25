import itertools
import collections
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
from Martini_PolyPly.itp_tool.itp_I import *
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
from numpy import sqrt, pi, cos, sin, dot, cross, arccos, degrees
from numpy import float as nfl
from numpy.linalg import norm
#=======================================================================================================================================================================
#                                                         GROMACS Specific Definitions      
#=======================================================================================================================================================================

# For adding a new term simply add it to the dictionary centers and key indexs to the dictionaries below together with the name between [].

centers = {   'moleculetype': [],
              'atoms': [0,5], 
              ('bonds', 1) : [0,1],
              ('angles', 1) : [0,1,2],
              ('constraints', 1): [0,1],
              ('dihedrals', 3): [0,1,2,3],
              ('dihedrals', 1): [0,1,2,3],
              'exclusions': [0,1,2,3],
              'virtual_sitesn': [0,2]}
settings ={
              'moleculetype':[0,1],
              'atoms':[1,2,3,4,6,7],
              ('bonds', 1) :[2,3,4], 
              ('angles',1):[3,4,5],
              ('constraints', 1):[2,3],
              ('dihedrals', 1):[4,5,6,7],
              'exclusions':[],
              'virtual_sitesn':[1,3,4]}

function ={ 'bonds':2, 
             'angles':3,
             'constraints':2,
             'dihedrals':4,
             'virtual_sitesn':0}

block_bonds={'PS' :{ 'PEO': '1 8000',
                     'P3HT': '1 8000',
                      'PP' : '1 8000'},
              'PEO':{  'PS':'1 8000',
                     'P3HT':'1 8000',
                      'PP' :'1 8000'},
              'P3HT':{'PEO': '1 8000',
                      'P3HT': '1 8000',
                       'PS' : '1 8000'}}

term_names=['moleculetype','atoms', 'bonds', 'angles', 'dihedrals', 'constraints', 'exclusions', 'virtual_sitesn']

# We could store the different format as a subdictionary and select based on the relevant function number in 
# the itp file. This would require modifcation of the write_itp function. 
# Use cases are : - dihedrals, - exclusions, virtual-sides and 

format_outfile={
                'bonds': '{:<5d} {:<5d} {:<2s} {:<8s} {:<8s}{}', 
                'angles': '{:<5d} {:<5d} {:<5d} {:<2s} {:<8s} {:<8s}{}',
                'dihedrals': '{:<5d} {:<5d} {:<5d} {:<5d} {:<2s} {:<8s}{:<8s}{:<1s}{}', 
                'atoms': '{:<5d} {:<5s} {:<1s} {:<5s} {:<3s} {:<1d} {:<8s} {:<3s}{}',
                'constraints': '{:<5d} {:<5d} {:<2s}{:<8s}{}','[': '{:<1s}{:<10s}{:<1s}{}',
                'moleculetype':'{:<5s} {:<1s}{}', 
                'exclusions': '{:<5d} {:<5d} {:<5d} {:<5d}{}',
                'virtual_sitesn':'{:<5d} {:<1s} {:<5d} {:<5s} {:<5s} {}'}

#=======================================================================================================================================================================
#                                                                         Summary of Functions
#=======================================================================================================================================================================

def line_up(new_centers):
    return([sorted(new_centers)[x][1] for x in np.arange(0,len(new_centers))])   

def move(center, count, n_atoms, offset):
    return(int(center) + n_atoms * count + offset)

def term_topology(key, term):
    if all([key != item for item in ['atoms', 'moleculetype', 'exclusions']]):
       return(centers[(key, int(term[function[key]]))], settings[(key, int(term[function[key]]))])
    else:
       return(centers[key], settings[key])

def repeat_term(term, key, n_trans, n_atoms, offset):
     count = 0
     new_terms = []
     max_atom = n_atoms * n_trans + offset
     center_indices, setting_indices = term_topology(key, term)

     while count < n_trans: 
          new_term = []
          [ new_term.append([x ,term[x]]) for x in setting_indices] 
          [ new_term.append([x, move(term[x], count, n_atoms, offset)]) for x in center_indices ]
          new_term = line_up(new_term)
 
          if all([ int(new_term[x]) <= max_atom for x in center_indices]):
             new_terms.append(new_term)
          count = count + 1

     return(new_terms)

def repeat_section(section, key, n_trans, n_atoms, offset):
       new_section = []
       for term in section:
           new_terms = repeat_term(term, key, n_trans, n_atoms, offset)
           [new_section.append(new_term) for new_term in new_terms]
       new_section=sorted(new_section)
       return(new_section)

def read_itp(name):
    itp = collections.OrderedDict({'moleculetype':[], 'atoms':[], 'bonds':[], 'angles':[], 'dihedrals':[], 'constraints':[], 'virtual_sitesn':[]})
    with open(name) as f:
         lines = f.readlines()
         for line in lines:
             words = line.split()
             if not any([ word in ';, \n, \r\n' for word in words]):
                if any([ word in '[ [ ]' for word in words]):
                   key = words[1]
                else:
                   add =  itp[key] + [line.replace('\n', '').split()]
                   itp.update({key:add})
    out_itp = collections.OrderedDict({})
    [ out_itp.update(collections.OrderedDict({key: value})) for key, value in itp.items() if len(value) != 0 ]
    return(out_itp)
          
def write_itp(text, outname):
    out_file = open(outname, 'w')
    for key in ['moleculetype', 'atoms', 'bonds', 'angles', 'dihedrals', 'constraints', 'virtual_sitesn']:
        if key in text:
           out_file.write('{:<1s}{:^18s}{:>1s}{}'.format('[',key,']','\n'))
           for line in text[key]:
               print(line)
               line.append('\n')
               out_file.write(str(format_outfile[key]).format(*line))

def itp_tool(itpfiles, n_mon, nexcl, outname, name): 
    block_count = 0 
    new_itp =collections.OrderedDict({'moleculetype':[], 'atoms':[], 'bonds':[], 'angles':[], 'dihedrals':[], 'constraints':[], 'virtual_sitesn':[]} )
    offset = 0
    n_atoms=0
    mon_itp = read_itp(itpfiles[0])
    nexcl = mon_itp["moleculetype"][0][1]
    new_itp.update({'moleculetype':[[name, nexcl]]})
  
    for name, n_trans in zip(itpfiles, n_mon):
           mon_itp = read_itp(name)
           n_atoms = len(mon_itp["atoms"])
           for key, section in mon_itp.items():                           
               if key != 'moleculetype':
                  add = new_itp[key] + repeat_section(section, key, n_trans, n_atoms, offset)
                  new_itp.update(collections.OrderedDict({key: add}))
               print(offset)
           offset += n_atoms * n_trans            
    out_itp = collections.OrderedDict({})
    [ out_itp.update({key: value}) for key, value in new_itp.items() if len(value) != 0 ]
    write_itp(out_itp, outname)
    return(None)

