import numpy as np
from tqdm import tqdm
#import scipy.spatiala
from cellgrid import capped_distance_array
import polyply

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
          self.molecule_list = []
          self.atom_list = []

      def add_atom(self, mol_name, mol_index, mol_atom_index,position,atom_index, top):

          atom_type = top.molecules[mol_name].atoms[mol_atom_index-1].parameters[0]

          try:
             self.positions.append(position)
          except AttributeError:
             self.positions = np.append(self.positions,position).reshape(-1,3)  

          self.atom_info.append((mol_name, mol_index, atom_type, mol_atom_index,atom_index))


      def remove_overlap(self,top, name, sol_names, cut_off, softness):
         print('probing for overlaps')
         if len(self.dist_matrix) == 0:
            self.distance_matrix(cut_off,top)

         bad_mol_indices = []
        
         for dist, info_A, info_B in self.dist_matrix:
          
             coefs, sigma, pot = polyply.structure_tool.force_field_tools.lookup_interaction_parameters(top, info_A[2], info_B[2], 'nonbond_params')
        
             if not polyply.structure_tool.force_field_tools.are_bonded_exception(info_A[3], info_B[3], info_A[0], top,'excl_list') and info_A[1] != info_B[1]:
               if info_A[0] != info_B[0] and dist < softness * sigma:
#                print(info_A[0],info_B[0])
                if info_A[0] != name:
                   bad_mol_indices.append(info_A[1])
                if info_B[0] != name:
                   bad_mol_indices.append(info_B[1])

         if len(bad_mol_indices) != 0:
            self.delete_molecules(bad_mol_indices) 

      def delete_molecules(self, bad_mol_indices):
          atom_indices = [ info[4] for info in self.atom_info if info[1] not in bad_mol_indices ]
#          print(atom_indices)
          self.positions = self.positions[atom_indices]
          new_indices = np.arange(0,len(atom_indices),1)
          self.atom_info[:] = [ (self.atom_info[i][0], self.atom_info[i][1], self.atom_info[i][2] , self.atom_info[i][3], j   ) for i,j in zip(atom_indices,new_indices)] 
          dist_indices = []
          count = 0

          for dist, info_A, info_B in self.dist_matrix:
              if info_A[1] not in bad_mol_indices or info_B[1] in bad_mol_indices:
                 dist_indices.append(count)
              count = count +1
          
          self.dist_matrix[:] = [ self.dist_matrix[i] for i in dist_indices] 

      @classmethod
      def from_gro_file(cls, name, top=None):
          lines = open(name).readlines()
          box = np.array([float(lines[-1].split()[0]), float(lines[-1].split()[1]), float(lines[-1].split()[2])])
          temp_traj = cls(box)
          big_res = False
          big_atoms = False    

          if top != None:
                 
              for molname, amount in top.composition:
                for i in np.arange(0,amount,1):
                    for atom in top.molecules[molname].atoms:
                        temp_traj.molecule_list.append((molname,i))
                        #print(int(atom.centers[0])-1)
                        temp_traj.atom_list.append(int(atom.centers[0]))
          atom_count = 0
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

             if big_res:
                res_count += 1
                res_index = res_count

             if res_index == 99999:
                big_res = True
                res_count = res_index

                  
             # if a topology is supplied, the types and resnames are matched to the topology
             if top != None:
                try:
                     mol_name   = temp_traj.molecule_list[atom_index-1][0]
                     mol_index  = temp_traj.molecule_list[atom_index-1][1]
                     mol_atom_index = temp_traj.atom_list[atom_index-1]
                     temp_traj.add_atom(mol_name, mol_index, mol_atom_index, point,atom_count,top)
 #                    print(atom_count)
                     atom_count += 1
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

      def distance_matrix(self,cut_off, top):
         
          traj_B=[]
          info_B=[]
          traj_B[:] = [ atom for atom in self.positions[:] ]
          info_B[:] = [ atom for atom in self.atom_info[:] ]

          self.dist_matrix=[]
          # This somewhat absurd loop structure is required to avoid calculating the distance pairs squared! 
          # One element is enough

          for atom_A, info_A in zip(self.positions,self.atom_info):
              traj_B[:] = [ atom for atom in traj_B[1:] ]
              info_B[:] = [ atom for atom in info_B[1:] ]
              ref_array = np.asarray(traj_B).reshape(-1,3) 
              if len(ref_array) - len(traj_B) != 0:
                 distances = capped_distance_array(atom_A.reshape(-1,3),ref_array, cut_off, self.box)
                 dist_list = [ (dist,info_A, info_B[key[1]]) for key, dist in zip(distances[0],distances[1]) if dist < cut_off and dist != 0.0] 
              else:
                 dist_list = []
              self.dist_matrix =  self.dist_matrix + dist_list      


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
