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
          
          def __init__(self, traj, cut_off):

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

          def add_atom(self, traj, new_atom, atom_type, molecule ,cut_off)
              flat_traj_A, info_A = traj.concatenate()
              ref_tree = scipy.spatial.ckdtree.cKDTree(new_atom.reshape(1,3))
              traj_tree = scipy.spatial.ckdtree.cKDTree(flat_traj_A)
              dist_mat = ref_tree.sparse_distance_matrix(traj_tree,cut_off)
              [ self.traj_dists.update({(atom_type, ):dist}) for key, dist in dist_mat.items()]


          def del_atom()



