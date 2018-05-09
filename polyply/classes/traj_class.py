import numpy as np
from tqdm import tqdm
#import scipy.spatiala
from cellgrid import capped_distance_array

class trajectory:
      '''
      An instance of the trajectory class stores the positions and topology information of the molecules 
      in the form of:
      '''

      def __init__(self, box=np.array([10.0,10.0,10.0])):
          self.box        = box
          self.positions  = []
          self.atom_info  = []
          self.n_atoms    = {}
          self.dist_matrix = []

      def add_atom(self, mol_name, mol_index, mol_atom_index,position,atom_index, top):
 #         print(mol_name)
          #print(top.molecules[mol_name].atoms[mol_atom_index-1].parameters[0])
          atom_type = top.molecules[mol_name].atoms[mol_atom_index-1].parameters[0]
          self.positions.append(position)
          self.atom_info.append((mol_name, mol_index, atom_type, mol_atom_index,atom_index))

      def remove_molecule(self, mol_name, index):
          try:
             for i, atom in enumerate(self.atoms):
                 if atom[1] == index:
                    del self.atoms[i]

          except (KeyError, IndexError):
             print("WARNING: Tried to delete molecule from trajectory but no molecule found.\n Will continue anyways.")

      @classmethod
      def from_gro_file(cls, name, top=None):
          lines = open(name).readlines()
          box = np.array([float(lines[-1].split()[0]), float(lines[-1].split()[1]), float(lines[-1].split()[2])])
          temp_traj = cls(box)
          big_res = False
          big_atoms = False    

          if top != None:
            molecule_list = [ ]
            atom_list = []
        
            for molname, amount in top.composition:
                for i in np.arange(0,amount,1):
                    for atom in top.molecules[molname].atoms:
                        molecule_list.append((molname,i))
                        #print(int(atom.centers[0])-1)
                        atom_list.append(int(atom.centers[0]))
        
          for line in lines[2:-1]:
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
             if top != None:
                try:
                     mol_name   = molecule_list[atom_index-1][0]
                     mol_index  = molecule_list[atom_index-1][1]
                     mol_atom_index = atom_list[atom_index-1]
                     temp_traj.add_atom(mol_name, mol_index, mol_atom_index, point,atom_index,top)
                except IndexError:
                     print("FATAL ERROR")
                     print("Check your topology and gro file format!")    
                     exit()
             # else we read the res_name as the molecule name 
             else:
                mol_name = "PSPEO"
                temp_traj.add_atom(mol_name, point, res_index-1, ff=None, atom_index=None)
          temp_traj.positions = np.asarray(temp_traj.positions).reshape(-1,3)
          return(temp_traj) 

      def distance_matrix(self,cut_off):
          traj_B=[]
          info_B=[]
          traj_B[:] = [ atom for atom in self.positions[:] ]
          info_B[:] = [ atom for atom in self.atom_info[:] ]

          self.dist_matrix=[]

          # This somewhat absurd loop structure is required to avoid calculating the distance pairs squared! 
          # One element is enough

          for atom_A, info_A in tqdm(zip(self.positions,self.atom_info)):
              traj_B[:] = [ atom for atom in traj_B[1:] ]
              info_B[:] = [ atom for atom in info_B[1:] ]
              ref_array = np.asarray(traj_B).reshape(-1,3)    
              distances = capped_distance_array(atom_A.reshape(-1,3),ref_array, cut_off, self.box)
              dist_list = [ (dist,info_A, info_B[key[1]]) for key, dist in zip(distances[0],distances[1]) if dist < cut_off and dist != 0.0] 
              self.dist_matrix =  self.dist_matrix + dist_list      

          #for item in self.dist_matrix:
          #    print(item)

      def add_atom_dist_matrix(self, new_atom, cut_off, mol_name, mol_index, atom_type, mol_atom_index):
          ref_tree = scipy.spatial.ckdtree.cKDTree(new_atom.reshape(1,3))
          traj_tree = scipy.spatial.ckdtree.cKDTree(self.positions)
          dist_mat = ref_tree.sparse_distance_matrix(traj_tree,cut_off)
          dist_list = [ (dist,(mol_name, mol_index, atom_type, mol_atom_index), self.atom_info[key[1]]) for key, dist in dist_mat.items() ]        

          self.dist_matrix =  self.dist_matrix + dist_list        


      def del_atom_dist_matrix(index):
              try:
                  del self.elemets[index]
              except IndexError:
                  print("WARNING: could not remove atom_pair!")
