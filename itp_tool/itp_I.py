import itertools
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

term_centers = {
              'moleculetype': [],
              'atoms': [0,2,5], 
              'bonds': [0,1],
              'angles': [0,1,2],
              'constraints': [0,1],
              'dihedrals': [0,1,2,3],
              'exclusions': [0,1,2,3],
              'virtual_sitesn': [0,2]}
term_keys ={
              'moleculetype':[0,1],
              'atoms':[1,3,4,6,7],
              'bonds':[2,3,4], 
              'angles':[3,4,5],
              'constraints':[2,3],
              'dihedrals':[4,5,6,7],
              'exclusions':[],
              'virtual_sitesn':[1,3,4]}

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
                'atoms': '{:<5d} {:<5s} {:<1d} {:<5s} {:<3s} {:<1d} {:<8s} {:<3s}{}',
                'constraints': '{:<5d} {:<5d} {:<2s}{:<8s}{}','[': '{:<1s}{:<10s}{:<1s}{}',
                'moleculetype':'{:<5s} {:<1s}{}', 
                'exclusions': '{:<5d} {:<5d} {:<5d} {:<5d}{}',
                'virtual_sitesn':'{:<5d} {:<1s} {:<5d} {:<5s} {:<5s} {}'}

#=======================================================================================================================================================================
#                                                                         Summary of Functions
#=======================================================================================================================================================================

def line_up(new_centers):
    return([sorted(new_centers)[x][1] for x in np.arange(0,len(new_centers))])   

def move(center, count, n_atoms):
    offset =0
    return(int(center) + n_atoms * count + offset)

def repeat_term(term, n_trans, n_atoms):
     count = 0
     new_terms = []
     max_atom = n_atoms * n_trans
     centers = term_centers[term[-1]]
     if term[-1] in '[ moleculetype ]':
        n_trans = 1
     while count < n_trans: 
          new_centers = [] 
          [new_centers.append([x, term[x]]) for x in term_keys[term[-1]]]
          [new_centers.append([x, move(term[x], count, n_atoms)]) for x in centers]
          new_term = line_up(new_centers)
          if all([ int(new_term[x]) <= max_atom for x in centers]):
             new_terms.append(new_term)
          count = count + 1
     return(new_terms)

def repeat_section(section, n_trans, n_atoms):
       new_section = []
       for term in section:
           new_terms = repeat_term(term, n_trans, n_atoms)
           [new_section.append(new_term) for new_term in new_terms]
       new_section=sorted(new_section)
       new_section.insert(0,[str(section[0][-1])])
       return(new_section)

def is_term(name, IDs):
    if name in IDs:
       return True
    else:
       return False

def get_sections_itp(itp, IDs):
     count=0
     all_terms=[]
     term_ID_list=[]
     n_IDs=len(IDs)
     while count < n_IDs:
       with open(itp) as f:
            lines=f.readlines()
            section=[]
            flag=False
            for line in lines:
                if not(any(is_term(word, [';', '\n', '\r\n']) for word in line.split())):
                  if len(line.split()) != 0:
                   if any(is_term(word, IDs) for word in line.split()):
                      flag = not flag
                      if len(section) == 0:
                         term_ID=(line.split()[1])
                         term_ID_list.append(line.split()[1])
                      else:
                         if term_ID in ['atoms']:
                            atoms=len(section)
                         all_terms.append(section)
                         section=[]          
                         break
                   elif flag:
                      term=line.replace('\n', '').split()
                      term.append(term_ID)
                      section.append(term)
            if len(section) != 0:
               all_terms.append(section)
       count = count + 1
       if term_ID in IDs:
          IDs.remove(term_ID)
     return(all_terms, atoms)

def write_itp(text, outname):
    out_file = open(outname, 'w')
    for item in text:
      for line in item:
        if str(line[0]) in '[ atoms, bonds, angles, dihedrals, constraints, moleculetype, exclusions, virtual_sitesn]':
           out_file.write('{:<1s}{:^18s}{:>1s}{}'.format('[',line[0],']','\n'))
           ID = line[0]
        else:
           line.append('\n')
       #    print(line)
           out_file.write(str(format_outfile[ID]).format(*line))

def itp_tool(itpfiles, n_mon, nexcl, outname): 
    block_count = 0 
    if len(itpfiles) > 1:
       new_itp=[]
       for name, n_trans in zip(itpfiles, n_mon):
           mon_itp, n_atoms = get_sections_itp(name, term_names)
           block = [ repeat_section(section, n_trans, n_atoms) for section in mon_itp ]
           [ new_itp.append(line) for line in block[2::len(block)]]
           block_count = n_trans * n_atoms
           print(block[0][1][0])
           new_itp = new_itp + block_bonds[block[0][1][0]]
       new_itp = np.insert(new_itp, ['[ moleculetype ]'])
       new_itp = np.insert(new_itp, [ nexcl + '3'])
    else:
       for name, n_trans in zip(itpfiles, n_mon):
           mon_itp, n_atoms = get_sections_itp(name, term_names)
           new_itp = [ repeat_section(section, n_trans, n_atoms) for section in mon_itp ]
    write_itp(new_itp, outname)
    return(None)

