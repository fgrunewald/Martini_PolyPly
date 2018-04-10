###########################################################################################
#                                                                                         #
#                            classes related to positioins                                #
#                                                                                         #
###########################################################################################

class trajectory:
      '''
      An instance of the trajectory class stores the positions and topology information of the molecules 
      in the form of:
      '''

      def __init__(self, box=np.array([10.0,10.0,10.0])):
          self.box       = box
          self.mol_pos   = {}
          self.mol_types = {}
          self.n_atoms   = 0
          self.n_mols    = {}

      def add_molecule(self, molname, positions, ff):
          types = []

          for index, position in enumerate(positions):
              atom_type = ff[molname]['atoms'][index]['typ']
              types += [atom_type]          
  
          try:
              self.mols[molname] += [positions]
              self.atom_types[molname] += [types]
       
          except KeyError:
              self.mol_pos[molname] = [positions]                         
              self.atom_types[molname] += [types]

          self.n_atoms = count_atoms(self.n_atoms, self.mols)
          self.n_mols  = count_mols (self.n_mols, self.mols)

     def remove_molecule(self, molname, index):
          try:
             self.n_atoms = self.n_atoms - len(self.mol_pos[molname][index])
             n = self.n_mols[molname] - 1     
             self.n_mols.update({molname:n})      

             del self.mol_pos[molname][index]
             del self.mol_types[molname][index]

          except KeyError, IndexError:
              print("WARNING: Tried to delete molecule from trajectory but no molecule found.\n Will continue anyways.")

     def count_atoms(init, mols):
         for name, molecules in mols.items():
             for infos in molecules:
                 init = init + len(infos[1]) 
         return(init)

     def count_mols(init, mols)
         for name, molecules in mols.items():
             init.update(len(molecules))             
         return(init)

     def concatenate(self):
         positions = np.concatenate([ pos  for molecule, pos_mols in self.mol_pos.items() for pos in pos_mols] 
         mol_types = [ (molecule, atom) for molecule, atoms in self.mol_types.items() for atom in atoms ]
           

class distance_matrix:

# Store as list???
# How do I determine if they are bonded??
          
          def __init__(self, traj, cut_off, topology):

              self.traj_dists={}

              flat_traj_A, info_A = traj.concatenate()
              flat_traj_B, info_B = traj.concatenate()

              for atom_A, info in zip(flat_traj_A[:-1], info_A):
           
                  flat_traj_B = np.delete(flat_traj_B,0,axis=0)
                  del info_B[0]
           
                  traj_tree = scipy.spatial.ckdtree.cKDTree(flat_traj_B)
                  ref_tree = scipy.spatial.ckdtree.cKDTree(atom_A.reshape(1,3))
                  dist_mat = ref_tree.sparse_distance_matrix(traj_tree,cut_off)
                  [ self.traj_dists.update({(info[0], mol_ids_A[count], atom_count, traj_info_B[key[1]][0], mol_ids_B[key[1]] , traj_info_B[key[1]][2]):dist}) for key, dist in dist_mat.items()] 

          def add_atom(self, traj, new_atom, atom_type, molecule ,cut_off):
              flat_traj_A, info_A = traj.concatenate()
              ref_tree = scipy.spatial.ckdtree.cKDTree(new_atom.reshape(1,3))
              traj_tree = scipy.spatial.ckdtree.cKDTree(flat_traj_A)
              dist_mat = ref_tree.sparse_distance_matrix(traj_tree,cut_off)
              [ self.traj_dists.update({(atom_type, ):dist}) for key, dist in dist_mat.items()]


          def del_atom(self,index)
              try:
                  del self[index]
              except IndexError:
                  print("WARNING: could not remove atom_pair!")

###########################################################################################################
#                                                                                                         #
#                                 Classes related to topology                                             #
#                                                                                                         #
###########################################################################################################

def strip_comments(line, topol_type):
    line = line.replace('\n', ' ').split()
    clean_line = []
    count = 0
    word = 'random'
    while True:
       if count < len(line):
          word = line[count]
          if word not in '[ ; ]':
             word = line[count]
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
      def __init__(self, key):
            self.key        = key
            self.centers = {}
            self.params = {}
            self.potential = {}
            self.term_types = []

      def add_term_type(self, term_type, centers, params, potential):
            self.centers.update({term_type:centers})
            self.params.update({term_type:params})
            self.term_types.append(term_type)
            self.potential.update({term_type:potential})

class topology_format:
      '''    
      This class can be used to store the entire topology informations of an MD program.
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
                      function = int(strip_comments(line)[1])

                   elif any([ word in '@params' for word in line.split()]):
                      params = change_type(strip_comments(line)[1:], int)

                   elif any([ word in '@centers' for word in line.split()]):
                      centers = change_type(strip_comments(line)[1:], int)

                   elif any([ word in '@potential' for word in line.split()]):
                      potential = strip_comments(line)[1:]

                   elif any([ word in '@add' for word in line.split()]):
                      term.add_term_type(function,centers, params, potential)
                      self.format.update({term.key:term})

      def defs(self, name):
          return(self.format[name].key, self.format[name].term_types, self.format[name].centers, self.format[name].params)



#########################################
#    Storing the topology informations
#########################################
                                       
class term_instance:
      
      def __init__(self, line, term_name, term_format):
      #    self.term_type     = term_type
          self.centers  = centers
          self.function_type = function
          self.paramters     = params
          self.potential     = potential

      def move(self, increment):
          new_centers = [atom + increment for atom in self.centers]
          new_term = term_instance(self.term_type)
          return(new_term.instance_direct(self.function_type, new_centers, self.paramters))

      @classmethod
      def from_line(cls, line, term_name, term_format, term_type):
         # self.term_type = term_type
          self.centers = line[term_format.centers]
          self.function_type = line[term_format.funct]
          self.params = line[term_format.params]
          self.potential = term_format.potential[self.function_type]


'''
These classes are initialized with specific class variables depeding on the MD program.
'''

class molecule:

      '''
       For the different MD codes we use alternative constructors in form of 
       of class methods. 
      '''
      def __init__(self, name):
          # a molecule needs to have at least a list of atoms and bonds and a name
          self.atoms = []
          self.bonds = []
          self.name = name

      def bonded_exclusions(self):
          self.mol_graph =  construct_mol_graph(self.bonds)
          exclusinos = {}
          for atom in self.atoms:
              exclusions.update({atom.centers: neighborhood(self.mol_graph, atom['n'], params['nexcl'])})
        
          return(exclusions) 
      
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
      def from_gromacs_lines(cls, lines, topology_format):
          self.name = line[0][1]
          indices =  [ (index, words[1]) for words in enumerate(lines) if line[1] in keywords ]
          
          for index, i in enumerate(indices):
              self.index[1] = []
              for term in line[index[0]:indices[i+1]]:
                  self.index[1] += term_instance.from_line(line, index[1], term_format)      

          self.exclusions = bonded_exclusions()

      @classmethod
      def from_namd_lines(cls):
          pass


class topology:

      def __init__(self, topol_format):
          
          self.composition = []
          self.molecules = []
          self.nonbondparams = []

      @classmethod
      def from_gromacs_topfile(cls, topol_file_name):
                      
          itp_files, toplines = find_itp_files(topol_file_name)
          lines = []
   
          for file_name in itp_files:
              lines += [ strip_comments(line, topol_type) for line in open(filename)]
              lines += ['EOF']
          
          

          def find_itp_files(topol_file_name):
              lines = open(topol_file_name)
              [line.replace('\n', ' ').replace('\"', ' ').split()[1]  for line in lines if any([ word in '#include' for word in line.split()])]
              return(file_list, lines)
  
          def split_lines(lines):
              return([(index, section) for index, word in strip_comments(lines.split()if any([ word in '#include' for word in line.split()]))  ]
          
          def extract_pairtypes(lines):
              pairtypes = {}
              for line in lines:
                  atom1, atom2, coef_I, coef_II = line.replace('\n',' ').split()
                  pairtypes.update{(atom1, atom2):{'coef_a','coef_b':}}   
              return(pairtypes)           

          def extract_nonbond(lines):
              nonbond_params = {}
              for line in lines:          
                  atom1, atom2, f, sigma, epsilon = line.replace('\n', '').split()
                  nonbond_params.update({(atom1, atom2): {'sigma':nfl(sigma), 'epsilon':nfl(epsilon)}})
              return(nonbond_params)     

          def extract_molecules(lines):
