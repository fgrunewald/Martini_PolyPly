import itertools
import collections
import numpy as np
#=======================================================================================================================================================================
#                                                         GROMACS Specific Definitions      
#=======================================================================================================================================================================

# For adding a new term simply add it to the dictionary centers and key indexs to the dictionaries below together with the name between [].

centers = {   'moleculetype': [],
              'atoms': [0,2,5],
              ('bonds', 1) : [0,1],
              ('bonds', 6) : [0,1],
              ('bonds', 2) : [0,1],
              ('position_restraints',1):[0],
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
              ('virtual_sitesn',2): [":"],
              ('virtual_sites3',4):[0,1,2,3],
              ('virtual_sites3',1):[0,1,2,3]}
settings ={
              'moleculetype':[0,1],
              'atoms':[1,3,4,6,7],
              ('bonds', 1) :[2,3,4],
              ('bonds', 2) : [2,3,4],
              ('bonds', 6) : [2,3,4],
              ('position_restraints',1):[1,2,3,4],
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
              ('virtual_sitesn',2):[],
              ('virtual_sites3',4):[4,5,6,7],
              ('virtual_sites3',1):[4,5,6]}

function ={ 'bonds':2,
             'angles':3,
             'constraints':2,
             'dihedrals':4,
             'pairs':2,
             'virtual_sitesn':1,
             'virtual_sites3':4,
             'position_restraints':1}


block_bonds={'PS' :{ 'PEO': '1 8000',
                     'P3HT': '1 8000',
                      'PP' : '1 8000'},
              'PEO':{  'PS':'1 8000',
                     'P3HT':'1 8000',
                      'PP' :'1 8000'},
              'P3HT':{'PEO': '1 8000',
                      'P3HT': '1 8000',
                       'PS' : '1 8000'}}

term_names=['moleculetype','atoms', 'bonds','position_restraints' ,'angles', 'dihedrals', 'constraints', 'exclusions', 'virtual_sitesn','virtual_sites3']

# We could store the different format as a subdictionary and select based on the relevant function number in 
# the itp file. This would require modifcation of the write_itp function. 
# Use cases are : - dihedrals, - exclusions, virtual-sides and 
format_outfile={
                ('bonds',1): '{:<5d} {:<5d} {:<2s} {:<8s} {:<8s}{}',
                ('angles',1): '{:<5d} {:<5d} {:<5d} {:<2s} {:<8s} {:<8s}{}',
                ('bonds',2): '{:<5d} {:<5d} {:<2s} {:<8s} {:<8s}{}',
                ('bonds',6): '{:<5d} {:<5d} {:<2s} {:<8s} {:<8s}{}',
                ('position_restraints',1): '{:<5d} {:<2s} {:<8s} {:<8s} {:<8s}{}',
                ('angles',2): '{:<5d} {:<5d} {:<5d} {:<2s} {:<8s} {:<8s}{}',
                ('angles',10): '{:<5d} {:<5d} {:<5d} {:<2s} {:<8s} {:<8s}{}',
                ('dihedrals',1): '{:<5d} {:<5d} {:<5d} {:<5d} {:<2s} {:<10s} {:<10s} {:<10s}{}',
                ('dihedrals',2):'{:<5d} {:<5d} {:<5d} {:<5d} {:<2s} {:<10s} {:<10s} {}',
                ('dihedrals',11):'{:<5d} {:<5d} {:<5d} {:<5d} {:<2s} {:<10s} {:<10s}  {:<10s} {:<10s} {:<10s} {:<10s} {}',
                ('dihedrals',9): '{:<5d} {:<5d} {:<5d} {:<5d} {:<2s} {:<10s} {:<10s} {:<10s}{}',
                'atoms': '{:<5d} {:<5s} {:<5d} {:<5s} {:<3s} {:<1d} {:<8s} {:<5s} {}',
                ('constraints',1): '{:<5d} {:<5d} {:<2s}{:<8s}{}','[': '{:<1s}{:<10s}{:<1s}{}',
                'moleculetype':'{:<5s} {:<1s}{}',
                'exclusions': '{:<5d} {:<5d} {}',
                ('virtual_sitesn',2):'{:<5d} {:<1s} {:<5d} {:<5d} {:<5d} {}',
                ('virtual_sites3',1):'{:<5d} {:<5d} {:<5d} {:<5d} {:<1s} {:<10s} {:<10s}{}',
                ('virtual_sites3',4):'{:<5d} {:<5d} {:<5d} {:<5d} {:<1s} {:<10s} {:<10s} {:<10s}{}',
                ('pairs',1):'{:<5d} {:<5d} {:<2s} {}'}

#=======================================================================================================================================================================
#                                                                         Summary of Functions
#=======================================================================================================================================================================

def line_up(new_centers):
    return([sorted(new_centers)[x][1] for x in np.arange(0,len(new_centers))])   

def move(center, count, n_atoms, offset):
    return(int(center) + n_atoms * count + offset)

def term_topology(key, term):
    # This takes care of define statements
    #print(term[0][0].split())
    if "#" in term[0][0].split():
        return (term, "define")
 
    if all([key != item for item in ['atoms', 'moleculetype', 'exclusions', 'virtual_sitesn']]):
       return(centers[(key, int(term[function[key]]))], settings[(key, int(term[function[key]]))])
    elif key == 'virtual_sitesn':
         #print(key)
         fidx = int(function[key])
         idxs = np.arange(0,len(term),1)
         cent_idxs = idxs[idxs != fidx]
         #print(fidx)
         #print(cent_idxs)
         return (cent_idxs, [fidx])
    else:
       return(centers[key], settings[key])

def repeat_term(term, key, n_trans, n_atoms, offset, max_atom,res_offset):
     count = 0
     new_terms = []
     center_indices, setting_indices = term_topology(key, term)
     #print(max_atom)
     #print(n_atoms)
     if setting_indices == "define":
        #print(center_indices) 
        return [[-1, center_indices]]

     while count < n_trans: 
          try:
              new_term = []
              [ new_term.append([x ,term[x]]) for x in setting_indices]  
              [ new_term.append([x, move(term[x], count, n_atoms, offset)]) for x in center_indices ]
              #print(new_term)
              #exit()
              new_term = line_up(new_term)
              #print(new_term)
          except IndexError:
              print("\n+++++++++++++++++++ Fatal Error +++++++++++++++++++++++")
              print("Check your itp file!")
              print("Too many or few parameters on the following line:")  
              print(term, '\n')
              exit()
          if key in '[ atoms ]':
             if new_term[center_indices[1]] > max_atom:
                #print(new_term) 
                print("\n++++++++++++++++ Fatal Error ++++++++++++++++++++++++")
                print("The largest charge group index exceeds the number")
                print("of atoms in the repeat unit. You cannot have more")
                print("charge groups than atoms. Check your input!\n")             
                exit()
          if all([ int(new_term[x]) <= max_atom for x in center_indices]):
             new_terms.append(new_term)
  #        print(new_term)
          count = count + 1

     # correction for couning the resids and the charge number
     if key == 'atoms':
        #print("term", " n_aotms*I "," offset "," i ")
        for i, term in enumerate(new_terms):
            #print(term[2],n_atoms*i ,res_offset, i)
            term[2] = term[2] - n_atoms*i - offset + i + res_offset
            #term[5] = term[5] - n_atoms*i - offset + i + 1

     return(new_terms)


# We need a special sorting algorithm to not sort around #defs

def check_interval(ndx,ifdef,endif):
    
    for idx, jdx in zip(ifdef,endif):
        if all([ndx >= idx  , ndx <= jdx]):
           return True
    return False   
    

def sort_section(section):
    
    sorted_section=[]
    ifdef=[]
    endif=[]
    
    if len(section) != 0:
       for i, term in enumerate(section):
           #print(term)
           try:
               #print(term[1][0])
               if  term[1][0] in ["#ifdef","#ifndef"] :
                   ifdef.append(i)
                 #  if term[1][0]  == "ifndef":
                 #     print("go here")
               elif "#endif" == term[1][0]:
                   endif.append(i)
            
           except TypeError:
               continue 
    #print(ifdef,endif)


    if len(ifdef) == 0 and len(section) != 0:
       #print("go") 
       #print(sorted(section)) 
       return sorted(section)

    else:
       #print("go here") 

       temp_sorted=[]
       for idx, term in enumerate(section):
           if not check_interval(idx, ifdef, endif):
              temp_sorted.append(term)
       sorted_section += sorted(temp_sorted)
       
       temp_sorted=[]
       for idx,jdx in zip(ifdef,endif):
           #print(idx,jdx)
           #print(section[idx])
           temp_sorted.append(section[idx])
           #print("ifdef",section[idx+1:jdx])
           temp_sorted += sorted(section[idx+1:jdx])
           temp_sorted.append(section[jdx])
           sorted_section += temp_sorted
           #print(sorted)

       return sorted_section

def repeat_section(section, key, n_trans, n_atoms, offset, max_atoms, res_offset):
       new_section = []
       for term in section:
           new_terms = repeat_term(term, key, n_trans, n_atoms, offset, max_atoms,res_offset)
           [new_section.append(new_term) for new_term in new_terms]
       #print(new_section)
       new_section=sort_section(new_section)
       #print(new_section)
       return(new_section)

def read_itp(name):
    itp = collections.OrderedDict({'moleculetype':[], 'atoms':[], 'bonds':[], 'angles':[], 'dihedrals':[], 'constraints':[],'position_restraints':[] ,'pairs':[],'virtual_sites3':[] ,'virtual_sitesn':[], 'exclusions':[]})
    with open(name) as f:
         lines = f.readlines()
         for line in lines:
             #print(line)
             words = line.split()
             if len(words) != 0:
               if not words[0] in ';, \n, \r\n':   
             #if not any([ word in ';, \n, \r\n' for word in words]):
                try:
                  if any([ word in '[ [ ]' for word in words]):
                     key = words[1]
                     #print(key)

                  elif key != 'exclusions':
                     add =  itp[key] + [line.replace('\n', '').split()]
                     itp.update({key:add})

                  elif key == 'exclusions':
                     sline = line.replace('\n', '').split()
                     #print(line)
                     for atom in sline[1:]:
                         #print(sline[0])
                         new_line = sline[0] + " " + atom
                         #print(new_line)
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
    #print(out_itp['bonds'])
    return(out_itp)

       
def write_itp(text, outname):
    out_file = open(outname, 'w')
    #print(text['bonds'])
    for key in ['moleculetype', 'atoms']:
      if key in text:
        out_file.write('{:<1s}{:^18s}{:>1s}{}'.format('[',key,']','\n'))
        for line in text[key]:
            line.append('\n')
            out_file.write(str(format_outfile[key]).format(*line))

    for key in ['bonds', 'angles', 'dihedrals', 'constraints','pairs','virtual_sites3','position_restraints']:
        if key in text:
           out_file.write('{:<1s}{:^22s}{:>1s}{}'.format('[',key,']','\n'))
           for line in text[key]:
               #print(line)
               if line[0] == -1:
                  #print(key) 
                  out_file.write(" ".join(line[1])+" \n")
               else: 
                  line.append('\n')
                  out_file.write(str(format_outfile[(key, int(line[function[key]]))]).format(*line))
    
    for key in ['exclusions', 'virtual_sitesn']:
        if key in text:
            out_file.write('{:<1s}{:^18s}{:>1s}{}'.format('[', key ,']','\n'))
            for line in text[key]:
                if line[0] == -1:
                   out_file.write(" ".join(line[1])+" \n")
                else: 
                    line.append('\n')
                    line = [ str(e) for e in line ]
                    out_file.write(" ".join(line))


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
    if len(itp['atoms']) != 0:
       last_res = itp['atoms'][-1][2] 
    for key, section in group.items():
        for term in section:
          #print(term)
          #print(offset)
          new_term = []
          center_indices, setting_indices = term_topology(key, term)
          [ new_term.append([x ,term[x]]) for x in setting_indices]
          if offset != 0:
             offset_new = -len(center_indices) + offset
             #print(len(center_indices))
             if key == 'atoms':
                offset_new = offset_new + 2
          else:
             offset_new = offset
          [ new_term.append([x, move(term[x], 0, 0, offset_new)]) for x in center_indices ]
          new_term = line_up(new_term)

          if key == 'atoms' and offset != 0:
             new_term[2] = last_res +1

          new_section = itp[key]
          new_section.append(new_term)
          itp.update({key: new_section})
    return(itp)    

def itp_tool(itpfiles, linkfile, n_mon, outname, name, term_info): 
    block_count = 0 
    new_itp =collections.OrderedDict({'moleculetype':[], 'atoms':[], 'bonds':[], 'angles':[], 'dihedrals':[], 'constraints':[],'position_restraints':[] ,'pairs':[] ,'virtual_sites3':[],'virtual_sitesn':[], 'exclusions':[]} )
    offset = 0
    n_atoms=0

    mon_itp = read_itp(itpfiles[0])
    nexcl = mon_itp["moleculetype"][0][1]
    new_itp.update({'moleculetype':[[name, nexcl]]})
    max_atoms = []
 
    for name, n_trans in zip(itpfiles, n_mon):
        mon_itp = read_itp(name)
        n_atoms = len(mon_itp["atoms"])
        try:
           max_atoms.append(n_atoms * n_trans + max_atoms[-1])
        except IndexError:
           max_atoms.append(n_atoms * n_trans)
        
    #print(n_atoms)
    if term_info != None:
       print("WARNING: The use of end-groups is to be deprecated.", '\n',
                        "Instead feed the end-group as monomer and", '\n',
                        "add corresponding link file.")
       new_itp = terminate(new_itp, term_info[0], 0)

#       if len(term_info) == 2:
#          atoms_last = len(read_itp(term_info[1])['atoms'])
#       else:
#          atoms_last = 0

       offset  = len(new_itp['atoms'])
#       max_atoms = [ n + len(new_itp['atoms']) + atoms_last for n in max_atoms ]  
       max_atoms = [ n + len(new_itp['atoms']) for n in max_atoms ]
   
    try:  
       res_offset = new_itp["atoms"][-1][2]
    except IndexError:
       res_offset = 0

    count=0
    for name, n_trans in zip(itpfiles, n_mon):
           mon_itp = read_itp(name)
           n_atoms = len(mon_itp["atoms"])
           for key, section in mon_itp.items():                           
               if key != 'moleculetype':
                  #print(max_atoms) 
                  add = new_itp[key] + repeat_section(section, key, n_trans, n_atoms, offset, max_atoms[count],res_offset)
                  new_itp.update(collections.OrderedDict({key: add}))
           #print(offset)
                  res_offset = new_itp["atoms"][-1][2]
           offset += n_atoms * n_trans            
           count = count + 1
    out_itp = collections.OrderedDict({})

    if linkfile != None:
       new_itp = add_links(new_itp, linkfile)

    if term_info != None and len(term_info) == 2:
       new_itp = terminate(new_itp, term_info[1], offset+1)
 
    [ out_itp.update({key: value}) for key, value in new_itp.items() if len(value) != 0 ]
    write_itp(out_itp, outname)
    return(None)

