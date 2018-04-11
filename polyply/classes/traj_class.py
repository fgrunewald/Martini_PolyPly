import numpy as np
from tqdm import tqdm

class trajectory:
      '''
      An instance of the trajectory class stores the positions and topology information of the molecules 
      in the form of:
   
      self.mol_pos contains a list of molecules each in turn with a numpy array of atoms
      self.mol_types contains a list of lists of atom-types, each list corresponds to one numpy array in self.mol_pos
      
      There are several ways to construct a trajectory:

      by adding molecules using add_molecule
      by adding atoms using add_atom
      by reading from a structure file in which case the molecules = residues unless topology information
      is provided 
      '''

      def __init__(self, box=np.array([10.0,10.0,10.0])):
          self.box       = box
          self.mol_pos   = {}
          self.atom_types = {}
          self.n_atoms   = {}
          self.n_mols    = {}

      def add_molecule(self, mol_name, positions=None, ff=None):
          types = []

          if len(positions) == 0:
             positions = np.array([])
    
          if ff != None:
             for index, position in enumerate(positions):
                 atom_type = ff[mol_name]['atoms'][index]['typ']
                 types += [atom_type]

          try:
              self.mol_pos[mol_name] += [positions.reshape(-1,3)]
              self.atom_types[mol_name] += [types]

          except KeyError:
              self.mol_pos.update({mol_name: [positions.reshape(-1,3)]})
              self.atom_types.update({mol_name: [types]})

          self.n_atoms = self.count_atoms(self.n_atoms, self.mol_pos)
          self.n_mols  = self.count_mols(self.n_mols, self.mol_pos)

      def add_atom(self, mol_name, positions, mol_index, ff=None, atom_index=None):
 
          if ff != None:
             atom_type = ff[mol_name]['atoms'][atom_index]['typ']
          else:
             atom_type = None

          try:
              self.mol_pos[mol_name][mol_index] = np.append(self.mol_pos[mol_name][mol_index], positions.reshape(-1,3), axis=0)
              self.atom_types[mol_name][mol_index] += [atom_type]

          except KeyError:
              self.add_molecule(mol_name, positions)
              
          except IndexError:
              self.mol_pos[mol_name].append(positions)
              self.atom_types[mol_name] += [[atom_type]]

          self.n_atoms = self.count_atoms(self.n_atoms, self.mol_pos)
          self.n_mols  = self.count_mols(self.n_mols, self.mol_pos)

      def remove_molecule(self, mol_name, index):
          try:
             self.n_atoms = self.n_atoms - len(self.mol_pos[mol_name][index])
             n = self.n_mols[mol_name] - 1
             self.n_mols.update({mol_name:n})

             del self.mol_pos[mol_name][index]
             del self.atom_types[mol_name][index]

          except (KeyError, IndexError):
             print("WARNING: Tried to delete molecule from trajectory but no molecule found.\n Will continue anyways.")

      def count_atoms(self, atom_counts, mol_pos):
          for mol_name, positions in mol_pos.items():
              atom_counts.update({mol_name:len(positions[0])})
          return(atom_counts)

      def count_mols(self, mol_counts, mol_pos):
          for mol_name, positions in mol_pos.items():
              mol_counts.update({mol_name:len(positions)})
          return(mol_counts)

 #     def concatenate(self):
 #        positions = np.concatenate([ pos  for molecule, pos_mols in self.mol_pos.items() for pos in pos_mols]
 #        mol_types = [ (molecule, atom) for molecule, atoms in self.mol_types.items() for atom in atoms ]

      @classmethod
      def from_gro_file(cls, name, ff=None):
          lines = open(name).readlines()

          temp_traj = cls(lines[-1])
          big_res = False
          big_atoms = False         
 
          for line in tqdm(lines[2:-1]):
             res_index = int(line.replace('\n', '')[0:5].strip())
             res_name = line.replace('\n', '')[5:10].strip()
             atom_name = line.replace('\n', '')[10:15].strip()
             atom_index = int(line.replace('\n', '')[15:20].strip())
             x = np.float(line.replace('\n', '')[20:28].strip())
             y = np.float(line.replace('\n', '')[28:36].strip())
             z = np.float(line.replace('\n', '')[36:44].strip())
             point = np.array([x,y,z])

             # This bit is to deal with gro files larger than 99999 atoms/residues
             if big_atoms:
                atom_count += 1
                atom_index = atom_count

             if big_res:
                res_count += 1
                res_index = res_count

             if res_index == 99999:
                big_res = True
                res_count = res_index

             if atom_index == 99999:
                big_atoms = True
                atom_count = atom_index
                  
             # if a topology is supplied, the types and resnames are matched to the topology
             if ff != None:
                print("Not implemtenred")
             # else we read the res_name as the molecule name 
             else:
                mol_name = res_name
                temp_traj.add_atom(mol_name, point, res_index-1, ff=None, atom_index=None)

          return(temp_traj)

