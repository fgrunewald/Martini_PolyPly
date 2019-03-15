import itertools
import collections
import numpy as np
#=======================================================================================================================================================================
#                                                         GROMACS Specific Definitions      
#=======================================================================================================================================================================

# For adding a new term simply add it to the dictionary centers and key indexs to the dictionaries below together with the name between [].

centers = {   'moleculetype': [],
              'atoms': [0,5], 
              ('bonds', 1) : [0,1],
              ('bonds', 2) : [0,1],
              ('angles', 1) : [0,1,2],
              ('angles', 2) : [0,1,2],
              ('angles', 10): [0,1,2],
              ('constraints', 1): [0,1],
              ('dihedrals', 9): [0,1,2,3],
              ('dihedrals', 11): [0,1,2,3],
              ('dihedrals', 1): [0,1,2,3],
              ('dihedrals', 3): [0,1,2,3],
              ('dihedrals', 2):[0,1,2,3],
              ('pairs',1):[0,1],
              'exclusions': [0, 1],
              ('virtual_sitesn',2): [0,2,3,4]}
settings ={
              'moleculetype':[0,1],
              'atoms':[1,2,3,4,6,7],
              ('bonds', 1) :[2,3,4],
              ('bonds', 2) : [2,3,4], 
              ('angles',1):[3,4,5],
              ('angles',2):[3,4,5],
              ('angles',10):[3, 4, 5],
              ('constraints', 1):[2,3],
              ('dihedrals', 1):[4,5,6,7],
              ('dihedrals',2):[4,5,6],
              ('dihedrals',11):[4,5,6,7,8,9,10],
              ('dihedrals', 3):[4,5,6,7,8,9,10],
              ('dihedrals',9):[4,5,6,7],
              ('pairs', 1):[2],
              'exclusions':[],
              ('virtual_sitesn',2):[1]}

function ={ 'bonds':2, 
             'angles':3,
             'constraints':2,
             'dihedrals':4,
             'pairs':2,
             'virtual_sitesn':1}

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
                ('bonds',1): '{:<5d} {:<5d} {:<2s} {:<8s} {:<8s}{}', 
                ('angles',1): '{:<5d} {:<5d} {:<5d} {:<2s} {:<8s} {:<8s}{}',
                ('bonds',2): '{:<5d} {:<5d} {:<2s} {:<8s} {:<8s}{}', 
                ('angles',2): '{:<5d} {:<5d} {:<5d} {:<2s} {:<8s} {:<8s}{}',
                ('angles',10): '{:<5d} {:<5d} {:<5d} {:<2s} {:<8s} {:<8s}{}',
                ('dihedrals',1): '{:<5d} {:<5d} {:<5d} {:<5d} {:<2s} {:<10s} {:<10s} {:<10s}{}', 
                ('dihedrals',2):'{:<5d} {:<5d} {:<5d} {:<5d} {:<2s} {:<10s} {:<10s} {}',
                ('dihedrals',11):'{:<5d} {:<5d} {:<5d} {:<5d} {:<2s} {:<10s} {:<10s}  {:<10s} {:<10s} {:<10s} {:<10s} {}',
                ('dihedrals',9): '{:<5d} {:<5d} {:<5d} {:<5d} {:<2s} {:<10s} {:<10s} {:<10s}{}', 
                'atoms': '{:<5d} {:<5s} {:<1s} {:<5s} {:<3s} {:<1d} {:<8s} {:<3s}{}',
                ('constraints',1): '{:<5d} {:<5d} {:<2s}{:<8s}{}','[': '{:<1s}{:<10s}{:<1s}{}',
                'moleculetype':'{:<5s} {:<1s}{}', 
                'exclusions': '{:<5d} {:<5d} {}',
                ('virtual_sitesn',2):'{:<5d} {:<1s} {:<5d} {:<5d} {:<5d} {}',
                ('pairs',1):'{:<5d} {:<5d} {:<2s} {}'}

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

def repeat_term(term, key, n_trans, n_atoms, offset, max_atom):
     count = 0
     new_terms = []
     center_indices, setting_indices = term_topology(key, term)
     #print(max_atom)
     #print(n_atoms)
     while count < n_trans: 
          try:
              new_term = []
              [ new_term.append([x ,term[x]]) for x in setting_indices] 
              [ new_term.append([x, move(term[x], count, n_atoms, offset)]) for x in center_indices ]
              new_term = line_up(new_term)
          except IndexError:
              print("\n+++++++++++++++++++ Fatal Error +++++++++++++++++++++++")
              print("Check your itp file!")
              print("Too many or few parameters on the following line:")  
              print(term, '\n')
              exit()
          if key in '[ atoms ]':
             if new_term[center_indices[1]] > max_atom:
                print("\n++++++++++++++++ Fatal Error ++++++++++++++++++++++++")
                print("The largest charge group index exceeds the number")
                print("of atoms in the repeat unit. You cannot have more")
                print("charge groups than atoms. Check your input!\n")             
                exit()
          if all([ int(new_term[x]) <= max_atom for x in center_indices]):
             new_terms.append(new_term)
  #        print(new_term)
          count = count + 1

     return(new_terms)

def repeat_section(section, key, n_trans, n_atoms, offset, max_atoms):
       new_section = []
       for term in section:
           new_terms = repeat_term(term, key, n_trans, n_atoms, offset, max_atoms)
           [new_section.append(new_term) for new_term in new_terms]
       new_section=sorted(new_section)
       return(new_section)

def read_itp(name):
    itp = collections.OrderedDict({'moleculetype':[], 'atoms':[], 'bonds':[], 'angles':[], 'dihedrals':[], 'constraints':[], 'pairs':[], 'virtual_sitesn':[], 'exclusions':[]})
    with open(name) as f:
         lines = f.readlines()
         for line in lines:
             words = line.split()
             if len(words) != 0:
               if not words[0] in ';, \n, \r\n':   
             #if not any([ word in ';, \n, \r\n' for word in words]):
                try:
                  if any([ word in '[ [ ]' for word in words]):
                     key = words[1]
                     #print(key)
                  else:
                     add =  itp[key] + [line.replace('\n', '').split()]
                     itp.update({key:add})       
                except (UnboundLocalError):
                       print("+++++++++++++ Error when reading the itp file ++++++++++++++++")
                       print("Check your format.")
                       print("There was something wrong with the section header!\n")
                       print("Note that there has to be a space between the [ and the section name.")
                       exit()
                except (KeyError):
                       print("+++++++++++++ Error when reading the itp file ++++++++++++++++")
                       print("Check your format.")
                       print("Your section type is currently not implemented.")
                       exit()
    out_itp = collections.OrderedDict({})
    [ out_itp.update(collections.OrderedDict({key: value})) for key, value in itp.items() if len(value) != 0 ]
 #   print(out_itp['atoms'])
    return(out_itp)

       
def write_itp(text, outname):
    out_file = open(outname, 'w')
    for key in ['moleculetype', 'atoms']:
      if key in text:
        out_file.write('{:<1s}{:^18s}{:>1s}{}'.format('[',key,']','\n'))
        for line in text[key]:
            line.append('\n')
            out_file.write(str(format_outfile[key]).format(*line))

    for key in ['bonds', 'angles', 'dihedrals', 'constraints','pairs','virtual_sitesn']:
        if key in text:
           out_file.write('{:<1s}{:^18s}{:>1s}{}'.format('[',key,']','\n'))
           for line in text[key]:
               line.append('\n')
               out_file.write(str(format_outfile[(key, int(line[function[key]]))]).format(*line))
    
    if 'exclusions' in text:
        out_file.write('{:<1s}{:^18s}{:>1s}{}'.format('[', 'exclusions',']','\n'))
        for line in text['exclusions']:
            line.append('\n')
            out_file.write(str(format_outfile['exclusions']).format(*line))



# The sole purpose of this function is to convert
# the centers to ints. So this can for sure be 
# handled smater in some way. 

def add_links(itp, linkfile):
    linkers = read_itp(linkfile)
    for key, section in linkers.items():
       for term in section:
          new_term = []
          center_indices, setting_indices = term_topology(key, term)
          [ new_term.append([x ,term[x]]) for x in setting_indices]
          [ new_term.append([x, int(term[x])]) for x in center_indices ]
          new_term = line_up(new_term)
          new_section = itp[key]
          new_section.append(new_term)
          itp.update({key: new_section})
    return(itp)
 
# We use the offset to automatically manipulate
# the terms of the end-group. Something similar
# could probably be done for the linker, so
# that one only has a single function. 

def terminate(itp, end_group_file, offset):
    group = read_itp(end_group_file)
    for key, section in group.items():
       for term in section:
          new_term = []
          center_indices, setting_indices = term_topology(key, term)
          [ new_term.append([x ,term[x]]) for x in setting_indices]
          if offset != 0:
             offset_new = -len(center_indices) + offset
             print(len(center_indices))
             if key == 'atoms':
                offset_new = offset_new + 1
          else:
             offset_new = offset
          [ new_term.append([x, move(term[x], 0, 0, offset_new)]) for x in center_indices ]
          new_term = line_up(new_term)
          new_section = itp[key]
          new_section.append(new_term)
          itp.update({key: new_section})
    return(itp)    

def itp_tool(itpfiles, linkfile, n_mon, outname, name, term_info): 
    block_count = 0 
    new_itp =collections.OrderedDict({'moleculetype':[], 'atoms':[], 'bonds':[], 'angles':[], 'dihedrals':[], 'constraints':[],'pairs':[] ,'virtual_sitesn':[], 'exclusions':[]} )
    offset = 0
    n_atoms=0

    mon_itp = read_itp(itpfiles[0])
    nexcl = mon_itp["moleculetype"][0][1]
    new_itp.update({'moleculetype':[[name, nexcl]]})
    max_atoms = []
 
    for name, n_trans in zip(itpfiles, n_mon):
        mon_itp = read_itp(name)
        n_atoms = len(mon_itp["atoms"])
        max_atoms.append(n_atoms * n_trans + sum(max_atoms))
 #   print(max_atoms[0])
    if term_info != None:
       print(term_info)
       new_itp = terminate(new_itp, term_info[0], 0)
       offset  = len(new_itp['atoms'])
       max_atoms = [ n + len(new_itp['atoms']) for n in max_atoms ]  

    count=0
    for name, n_trans in zip(itpfiles, n_mon):
           mon_itp = read_itp(name)
           n_atoms = len(mon_itp["atoms"])
           for key, section in mon_itp.items():                           
               if key != 'moleculetype':
                  add = new_itp[key] + repeat_section(section, key, n_trans, n_atoms, offset, max_atoms[count])
                  new_itp.update(collections.OrderedDict({key: add}))
                #print(offset)
           offset += n_atoms * n_trans            
           count = count + 1
    out_itp = collections.OrderedDict({})

    if linkfile != None:
       new_itp = add_links(new_itp, linkfile)

#    try:
    new_itp = terminate(new_itp, term_info[1], offset+1)
 #   except (IndexError,TypeError):
   #     print('WARNING: No second terminal found. Are you sure you want to leave the polymer unterminated?')
 
    [ out_itp.update({key: value}) for key, value in new_itp.items() if len(value) != 0 ]
    write_itp(out_itp, outname)
    return(None)

