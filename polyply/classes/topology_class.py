###########################################################################################################
#                                                                                                         #
#                                 Classes related to topology                                             #
#                                                                                                         #
###########################################################################################################

import networkx as nx


#################
# File parseing
#################


def strip_comments(line):
     line = line.replace('\n', ' ').split()
     clean_line = []
     count = 0
     word = 'random'
     while True:
        if count < len(line):
           word = line[count]
           if word != ';':
              if word != ']' and word != '[':
                 clean_line += [word]
              count = count + 1
           else:
              break
        else:
           break
     return(clean_line)


def change_type(line,py_type):
    return([ py_type(word) for word in line])


class read_top_file:
      '''
      This class stores potential input parsers for different MD programs. 
      At the moment only gromacs topologies can be parsed.
      '''
      @classmethod
      def gromacs(cls,topol_file_name):
     
         import os
         import re
     
         #############################################
         ### Start from backwards.py by T. Wassener
         #############################################
     
         includePattern = re.compile('#include "(.*)"')
     
         # Gromacs force field directory
         gmxlib = os.environ.get("GMXLIB")
         if not gmxlib:
               gmxdat = os.environ.get("GMXDATA")
               if gmxdat:
                   gmxlib = os.path.join(gmxdat,"gromacs","top")
               else:
                   gmxlib="."
         
         def reciter(filename):
           # Set the directory of the filename so we know where to expect #included files
           dir = os.path.dirname(filename)
         
           # Iterate over the lines
           lines = open(filename).readlines()   # we need to keep track on when a file ends
           lines += ['EOF']                     # this is a small modefication   
           for line in lines:
         
               # Check for an #include statement; yield the line if there is none
               if line.strip().startswith("#include"):
         
                   # Extract the #include filename             
                   matches = re.findall(includePattern,line)
         
                   if matches:
                       fr = matches[0]
         
                       if not os.path.exists(fr):
                           fr = os.path.join(dir,matches[0])
         
                       if not os.path.exists(fr):
                           fr = os.path.join(gmxlib,matches[0])
         
                       if not os.path.exists(fr):
                           yield "; " + line + " ; File not found\n"
                       else:
                           for j in reciter(fr):
                               yield j
               else:
                   yield strip_comments(line) # small modification to get rid of comments
         
                    
         #############################################
         ### End from backwards.py 
         #############################################

         
         lines = [ line for line in reciter(topol_file_name) if len(line) != 0]

         return lines

  
##################
#     FORMAT
##################      


class term_format:

      '''
      Term format is a general class for handeling the format of terms. 
      Each term can have serveral function attributes which are stored
      in this class as well. 
      '''

      def __init__(self, term_name, ftl=0):
           
            # the name of the term (e.g. dihedrals, angles, bonds)
            self.term_name  = term_name
           
            # field where the function type is located
            self.identifier = ftl
           
            # fields where the atom numbers are located
            self.centers = {}
           
            # fields where the parameters are located except 
            # function type 
            self.params = {}
           
            # the potential which applies to this function type
            self.potential = {}
           
            # the function types which are available 
            self.term_types = []

      def add_term_type(self, function, function_pos, potential, centers, params):
            self.identifier = function_pos
            self.centers.update({function:centers})
            self.params.update({function:params})
            self.term_types.append(function)
            self.potential.update({function:potential})

class topology_format:
      '''    
      This class can be used to store the entire topology informations of an MD program.
      It reads a simple text file with six commands defining the structure of the topology.

      The format file takes a ; as comment character. 
      
      You can use the following fields to feed the format class:
      term-name   function-location  function-type   potential     atom-centers    parameters  
      
      Each section can have serveral functions. 

      '''
      def __init__(self, name, filename):
          self.md_code = name
          self.format = {}
        
          with open(filename) as f:
             lines = f.readlines()
         
             for line in lines:
                 line = ''.join(strip_comments(line))
                 line = line.split('&')
         
                 if line != ['']:
                    print(line)
                    
                    try:
                        term = self.format[line[0]]
                    except KeyError:
                        term = term_format(line[0])
                    
                    function = line[2]
                    if line[1] != 'dum':
                       function_pos = int(line[1])
                    else:
                       function_pos = line[1]
                    potential = line[3]
                    centers = list(map(int,list(line[4])))
                    params  = list(map(int,list(line[5])))
                    term.add_term_type(function, function_pos, potential, centers, params)
                    self.format.update({line[0]:term})


      def defs(self, name):
          return(self.format[name].key, self.format[name].term_types, self.format[name].centers, self.format[name].params)

#########################################
#    Storing the topology informations
#########################################

class term_instance:
      '''
      These classes are initialized with specific class variables depeding on the MD program.
      This will be our eierlegende Woll-milch-sau (ELWMS) object
      '''

      def __init__(self, centers, term_type, params, potential):
      #    self.term_type     = term_type
          self.centers       = centers
          self.function_type = term_type
          self.parameters     = params
          self.potential     = potential

      def move(self, increment):
          new_centers = [atom + increment for atom in self.centers]
          new_term = term_instance(self.term_type)
          return(new_term.instance_direct(self.function_type, new_centers, self.paramters))

      @classmethod
      def from_line(cls, line, term_name, term_format):
          try:
             term_type = line[term_format.identifier]
          except TypeError:
             term_type = term_format.identifier

          #print(line)

          centers = tuple([line[int(element)] for element in term_format.centers[term_type]])
          params = [line[int(element)] for element in term_format.params[term_type]]
          potential = term_format.potential[term_type]    

          return(cls(centers, term_type, params, potential))

class top_base:

      def __init__(self):
          pass

      def get_indices(lines, keywords):
          return([(index, word) for index, line in enumerate(lines) for word in line if word in keywords ])


class molecule(top_base):

      '''
       For the different MD codes we use alternative constructors in form of 
       of class methods. 
      '''
      def __init__(self, name, excl):
          # a molecule needs to have at least a list of atoms and bonds and a name
          self.name = name
          self.excl = int(excl)

      def add_potential_term(self, name, line, term_format):
 #         print(name)
          try:
             getattr(self,name).append(term_instance.from_line(line, name, term_format))
          except AttributeError:
             setattr(self,name,[])
             getattr(self,name).append(term_instance.from_line(line, name, term_format))   

      def add_residue(self, resname, program='gromacs', resdir='def'):
          
          lines = read_itp(resdir, program)
          program_parser = 'from_' + program + 'lines'         
          residue = self.getattr(self,program_paser)(lines)
          
             


 
      def construct_mol_graph(self):
          G = nx.Graph()
          edges = [ (int(entry.centers[0]), int(entry.centers[1])) for entry in self.bonds]
          G.add_edges_from(edges)
          if not nx.is_connected(G): 
             print("You fed me an ill ordered topology: I'll proceed and die. x_x")
             exit()
          print(G.edges)
          print(G.nodes)
          mol_tree = nx.dfs_tree(G,1)
          print(mol_tree.edges())
          print("$$$$$$$$$$$$$")
          return(G)

      def neighborhood(self, node, n):
          G = self.mol_graph
          # Adobted from: https://stackoverflow.com/questions/22742754/finding-the-n-degree-neighborhood-of-a-node    
          path_lengths = nx.single_source_dijkstra_path_length(G, node)
          neighbours=[ node for node, length in path_lengths.items() if length <= n]
          return(neighbours)
 
      def bonded_exclusions(self):
          try:
             self.mol_graph =  self.construct_mol_graph()
             self.excl_list = {}
             self.excl_14_list = {}
       
             for atom in self.atoms:
                 self.excl_list.update({atom.centers[0]: self.neighborhood(int(atom.centers[0]), self.excl)})
                 self.excl_14_list.update({atom.centers[0]: self.neighborhood(int(atom.centers[0]),4)})
          except AttributeError:
             self.excl_list={}
             self.excl_14_list = {}
 
      def convert_constraints(self):
          new_bonds = []
          try:
             new_bonds = [ term_instance(term.centers, '1', [term.parameters[0], 9000], 'SimpleHarmonic')  for term in self.constraints ]
          except AttributeError:
             return
          try:
             self.bonds += new_bonds
             return
          except AttributeError:
             setattr(self,'bonds',new_bonds)
             return

      @classmethod
      def from_gromacs_lines(cls, lines, topology_format,end):
          keywords = ['atoms','bonds','angles','dihedrals','constraints','virtualsides'] + [end]
          mol = cls(lines[1][0], lines[1][1])
          indices =  cls.get_indices(lines,keywords)
          for i, index in enumerate(indices[:-1]):
              for line in lines[index[0]+1:indices[i+1][0]]:
                  mol.add_potential_term(index[1], line, topology_format.format[index[1]])
          mol.convert_constraints()
          mol.bonded_exclusions()
          return(mol)

      @classmethod
      def from_namd_lines(cls, lines, topology_format):
          pass

class topology(top_base):

      def __init__(self, name):

          self.name = name
          self.composition = []
          self.molecules = {}
          self.nonbondparams = {}
          self.nonbondparams_14 = {}
          self.bondtypes = {}
          self.defaults = {}

      def add_nonmol_params(self, lines, topology_format):
       #   print(lines)
      #    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
          term_name = lines[0][0]
          setattr(self,term_name,{})
          for line in lines[1:-1]:
              term = term_instance.from_line(line, term_name, topology_format.format[term_name])   
              getattr(self,term_name).update({term.centers:term})

      @classmethod
      def from_gromacs_topfile(cls, topol_file_name, topology_format):


          top_sections = 'EOF defaults molecules nonbond_params moleculetype bondtypes constrainttypes angletypes dihedraltypes pairtypes'.split()
         
          print(topol_file_name)
          lines = read_top_file.gromacs(topol_file_name)

          sys_index = topology.get_indices(lines, ['system'])[0][0]
          sys_name = lines[sys_index+1]
          top = topology(sys_name)
          indices = topology.get_indices(lines, top_sections)
          #print(lines)
          for i, index in enumerate(indices):
              # we now load the topology accoding to the section in the GROMACS format     
          
              if   index[1] == 'moleculetype':
                   mol = molecule.from_gromacs_lines(lines[index[0]:indices[i+1][0]+1],topology_format, indices[i+1][1])
                   top.molecules.update({mol.name:mol})

              elif index[1]  == 'molecules':
                   for line in lines[index[0]+1:indices[i+1][0]]:
                       top.composition.append((line[0],int(line[1])))

              elif index[1] == 'defaults':
                   top.defaults.update({'LJ':lines[index[0]+1][1],'COUL':lines[index[0]+1][0]})

              elif index[1] != 'EOF' and index[1] != 'system':
                   #print(lines[index[0]:indices[i+1][0]])
                   print(index,indices[i+1])
                   top.add_nonmol_params(lines[index[0]:indices[i+1][0]],topology_format)
             
              else:
                   pass
          return(top)
