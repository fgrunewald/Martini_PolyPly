###########################################################################################################
#                                                                                                         #
#                                 Classes related to topology                                             #
#                                                                                                         #
###########################################################################################################

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


##################
#     FORMAT
##################      


class term_format:
      '''
      Term format is a general class for handeling the format of terms. 
      '''
      def __init__(self, key, idf=0):
            self.key        = key
            self.identifier = idf
            self.centers = {}
            self.params = {}
            self.potential = {}
            self.term_types = []

      def add_term_type(self, term_type, funct_pos, centers, params, potential):
            self.identifier = funct_pos
            self.centers.update({term_type:centers})
            self.params.update({term_type:params})
            self.term_types.append(term_type)
            self.potential.update({term_type:potential})

class topology_format:
      '''    
      This class can be used to store the entire topology informations of an MD program.
      It reads a simple text file with six commands defining the structure of the topology.
 
      [ section ] indicates what we are dealing with (i.e. bonds, angles, dihedrals, nonbonded paramters etc.)

      @function  tells us on what field if at all a function type is specified
      @params    gives the fields assciated to information that are not atom centers
      @centers   gives the fields corresponding to atom centers
      @potential defines what the potential function is we are dealing with (e.g. harmonic, LJ, cosine, periodic etc. )
      @add indicates the end of a function section

      Each section can have serveral functions. 

      '''
      def __init__(self, name, filename):
          self.md_code = name
          self.format = {}

          with open(filename) as f:
               lines = f.readlines()
               term_name = 'start'
               function = None
               keyword = None
               for line in lines:
                   if any([ word in '[ [ ]' for word in line.split()]):
                      term_name = line.replace('\n', '').split()[1]
                      term = term_format(term_name)

                   elif any([ word in '@function' for word in line.split()]):
                      try:
                         function = int(strip_comments(line)[1])
                      except ValueError:
                         function = strip_comments(line)[1]

                   elif any([ word in '@params' for word in line.split()]):
                      params = change_type(strip_comments(line)[1:], int)

                   elif any([ word in '@centers' for word in line.split()]):
                      centers = change_type(strip_comments(line)[1:], int)

                   elif any([ word in '@potential' for word in line.split()]):
                      potential = strip_comments(line)[1]
                      pot_type = strip_comments(line)[2]

                   elif any([ word in '@add' for word in line.split()]):
                      term.add_term_type(pot_type, function, centers, params, potential)
                      self.format.update({term.key:term})

      def defs(self, name):
          return(self.format[name].key, self.format[name].term_types, self.format[name].centers, self.format[name].params)

#########################################
#    Storing the topology informations
#########################################

class term_instance:
      '''
      These classes are initialized with specific class variables depeding on the MD program.
      '''

      def __init__(self, centers, term_type, params, potential):
      #    self.term_type     = term_type
          self.centers       = centers
          self.function_type = term_type
          self.paramters     = params
          self.potential     = potential

      def move(self, increment):
          new_centers = [atom + increment for atom in self.centers]
          new_term = term_instance(self.term_type)
          return(new_term.instance_direct(self.function_type, new_centers, self.paramters))

      @classmethod
      def from_line(cls, line, term_name, term_format):
          try:
             print(term_format.identifier)
             print(line)
             term_type = line[term_format.identifier]
          except TypeError:
             term_type = term_format.identifier

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
      def __init__(self, name):
          # a molecule needs to have at least a list of atoms and bonds and a name
          self.atoms = []
          self.bonds = []
          self.name = name

      def add_potential_term(self, name, line, term_format):
          print(name)
          if name in locals():
             self.name.append(term_instance.from_line(line, name, term_format))
          else:
             self.name = []
             self.name.append(term_instance.from_line(line, name, term_format))   

      def bonded_exclusions(self):
          # We only construct these once
          if self.exclusions not in locals():
             self.mol_graph =  construct_mol_graph(self.bonds)
             self.exclusinos = {}
             for atom in self.atoms:
                 exclusions.update({atom.centers: neighborhood(self.mol_graph, atom['n'], params['nexcl'])})

          def construct_mol_graph(bonds):
              G = nx.Graph()
              edges = [ (entry['pairs'][0], entry['pairs'][1]) for entry in bonds]
              G.add_edges_from(edges)
              return(G)

          def neighborhood(G, node, n):
              # Adobted from: https://stackoverflow.com/questions/22742754/finding-the-n-degree-neighborhood-of-a-node    
              path_lengths = nx.single_source_dijkstra_path_length(G, node)
              neighbours=[ node for node, length in path_lengths.items() if length <= n]
              return(neighbours)

      @classmethod
      def from_gromacs_lines(cls, lines, topology_format,end):
          keywords = ['atoms','bonds','angles','dihedrals','constraints'] + [end]
          mol = cls(lines[1][1])
          indices =  cls.get_indices(lines,keywords)
          for i, index in enumerate(indices):
              for line in lines[index[0]+2:indices[i+1][0]]:
                  print(line)
                  mol.add_potential_term(index[1], line, topology_format.format[index[1]])

          mol.bonded_exclusions
          return(mol)

      @classmethod
      def from_namd_lines(cls, lines, topology_format):
          pass

class topology(top_base):

      def __init__(self, name):

          self.name = name
          self.composition = {}
          self.molecules = {}
          self.nonbondparams = {}
          self.nonbondparams_14 = {}
          self.bondtypes = {}
          self.defaults = {}

      def add_nonmol_params(self, lines, topology_format):
          #print(lines)
          term_name = lines[0][0]
          specification = {}
          for line in lines[1:-1]:
              print(line)
              term = term_instance.from_line(line, term_name, topology_format.format[term_name])   
             # print(term.centers)    
              specification.update({term.centers:term})
          self.term_name = specification

#      @staticmethod
#      def get_indices(lines, keywords):
#          return([(index, word) for index, line in enumerate(lines) for word in line if word in keywords ])

      @classmethod
      def from_gromacs_topfile(cls, topol_file_name, topology_format):

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

          top_sections = '[defaults, molecules, nonbond_params, moleculetype, bondtypes, constrainttypes, angletypes, dihedraltypes, pairtypes ]'
          lines = [ line for line in reciter(topol_file_name) if len(line) != 0]
          sys_index = topology.get_indices(lines, ['system'])[0][0]
          sys_name = lines[sys_index+1]
          top = topology(sys_name)
          indices = topology.get_indices(lines, top_sections)

          for i, index in enumerate(indices):
              # we now load the topology accoding to the section in the GROMACS format     
          
              if   index[1] == 'moleculetype':
                   mol = molecule.from_gromacs_lines(lines[index[0]:indices[i+1][0]+1],topology_format, indices[i+1][1])
                   top.molecules.update({mol.name:mol})

              elif index[1]  == 'molecules':
                   for line in lines[index:indices[j+1]]:
                        try:
                           n_mol = top.composition[line[0]] + int(line[1])
                           top.composition.update({line[0]:n_mol})
                        except KeyError:
                           top.composition.update({line[0]:int(line[1])})

              elif index[1] == 'defaults':
                   top.defaults.update({'LJ':lines[index[0]+1][1]})

              elif index[1] != 'EOF' and index[1] != 'system':
                   top.add_nonmol_params(lines[index[0]:indices[i+1][0]],topology_format)
             
              else:
                   pass
