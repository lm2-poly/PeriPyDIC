#-*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrick.diehl@polymtl.ca

import logging
import numpy as np
import util.neighbor
from scipy import linalg 
np.set_printoptions(threshold='nan')
import sys

logger = logging.getLogger(__name__)

class PD_problem():

    def __init__(self, deck):
        # Import initial data
       
        # NeighborSearch
        self.neighbors = util.neighbor.NeighborSearch(deck)

        # Compute the external force density applied on each node in a vector.
        self.compute_b(deck)

        # Compute the weighted volume for each node in a vector.
        self.weighted_function(deck)
        
        # Initialize the vector saving the data
        self.y = np.zeros((deck.num_nodes, deck.dim, deck.time_steps),dtype=np.float64)
        self.y[:,:,0] = deck.geometry.nodes[:,:]

        self.force_int = np.zeros((deck.num_nodes, deck.dim, deck.time_steps),dtype=np.float64)
        self.ext = np.zeros( ( deck.num_nodes, deck.num_nodes, deck.time_steps ),dtype=np.float64 )

    # Creates a loading vector b which describes the force or displacement applied on each node
    # at any time step
    def compute_b(self, deck):
        #Force density as a 3D array with the last dimension for the time steps
        self.b = np.zeros((deck.num_nodes, deck.dim, deck.time_steps),dtype=np.float64)
        for t_n in range(1, deck.time_steps):
            for con in deck.conditions:
                if con.type == "Force":
                    #Madenci approach
                    for i in con.id:
                        # x direction
                        if con.direction == 1:
                            self.b[int(i), 0, t_n] = self.shape_loading( deck, t_n , con , i )
                        # y direction
                        if con.direction == 2:
                            self.b[int(i),  1 , t_n] = self.shape_loading( deck, t_n , con , i )
                        # z direction
                        if con.direction == 3:
                            self.b[int(i), 2, t_n] = self.shape_loading( deck, t_n , con , i )
        #print self.b

    # Provide the loading shape to use to compute the loading vector b
    def shape_loading(self, deck, t_n, con, i):
        Time_t = deck.delta_t*(t_n)
        if deck.shape_type == "Ramp":
            force_density = con.value / deck.geometry.volumes[int(i)]
            if Time_t <= deck.shape_values[0]:
                result = (force_density*Time_t)/deck.shape_values[0]
                return result
            elif Time_t > deck.shape_values[0] and Time_t <= deck.shape_values[1]:
                result = force_density
                return result
            elif Time_t > deck.shape_values[1] and Time_t <= deck.shape_values[2]:
                result = force_density - force_density*(Time_t - deck.shape_values[1])/(deck.shape_values[2] - deck.shape_values[1])
                return result
            else:
                result = 0
                return result
        else:
            logger.error("Error in problem.py: Shape of BC unknown, please use Ramp.")

    # Computes the weights for each PD node
    def weighted_function(self, deck):
        self.weighted_volume = np.zeros((deck.num_nodes),dtype=np.float64)
        for i in range(0, deck.num_nodes):
            index_x_family = self.neighbors.get_index_x_family(i)
            for p in index_x_family:
                X = deck.geometry.nodes[p,:] - deck.geometry.nodes[i,:]
                self.weighted_volume[i] += deck.influence_function * (np.linalg.norm(X))**2 * deck.geometry.volumes[p]

    # Computes the internal force density
    def compute_f(self, y, deck, t_n):
        # Choice of the material class
        if deck.material_type == "Elastic":
            from materials.elastic import Elastic_material
            self.mat_class = Elastic_material( deck, self, y )
            self.update_force_data(self.mat_class, t_n)
            self.update_ext_state_data(self.mat_class, t_n)
            
        elif deck.material_type == "Viscoelastic":
            from materials.viscoelastic import Viscoelastic_material
            self.mat_class = Viscoelastic_material( deck, self, y, t_n)
            self.update_force_data(self.mat_class, t_n)
            self.update_ext_state_data(self.mat_class, t_n)
            self.update_ext_state_visco_data(self.mat_class, t_n)
        else:
            logger.error("Error in problem.py: Material type unknown, please use Elastic or Viscoelastic.")
        
        internal_force = self.mat_class.f_int
        #internal_force = np.reshape(internal_force, (deck.num_nodes * deck.dim,-1) ) 
        return internal_force
    
    # Computes the residual vector used in the quasi_static_solver function
    def compute_residual(self, y, deck, t_n):
        residual = np.zeros((deck.num_nodes, deck.dim),dtype=np.float64)
        internal_force = self.compute_f(y, deck, t_n)
        for con in deck.conditions:
            if con.type == "Displacement":
                for id_node in con.id:
                    # x direction
                    if con.direction == 1:
                        y[int(id_node),0] = deck.geometry.nodes[int(id_node),0] + con.value
                    # y direction
                    if con.direction == 2:
                        y[int(id_node),1] = deck.geometry.nodes[int(id_node),1] + con.value
                    # z direction
                    if con.direction == 3:
                        y[int(id_node),2] = deck.geometry.nodes[int(id_node),2] + con.value
        for i in range(0,deck.num_nodes):
            found = False
            for con in deck.conditions:
                if con.type == "Displacement":
                    if i in con.id:
                        found = True
            if found == False:  
                residual[i,:] = internal_force[i,:] + self.b[i,:, t_n]               
        return residual

    def compute_jacobian(self, y, deck, t_n, perturbation_factor):
        eps = perturbation_factor * deck.delta_X
        jacobian = np.zeros((deck.num_nodes * deck.dim , deck.num_nodes * deck.dim),dtype=np.float64)
            
        ids = []   
        for con in deck.conditions:
            if con.type == "Displacement":
                for i in con.id:
                    ids.append(i)

        for i in range(0, deck.num_nodes):
            traversal_list = np.append([i],self.neighbors.get_index_x_family(i))
            for j in traversal_list :
                for r in range(0, deck.dim):
                    eps_vector = np.zeros((deck.num_nodes , deck.dim),dtype=np.float64)
                    eps_vector[j,r] = eps
                    force_int_p = self.compute_f(y + eps_vector, deck, t_n)[i,:]
                    #force_int_p = np.array([[1],[2],[3],[4],[5],[6]])
                    #force_int_p = np.array([[1,10],[2,20],[3,30],[4,40],[5,50],[6,60]])
                    force_int_m = self.compute_f(y - eps_vector, deck, t_n)[i,:]
                    #force_int_m = np.array([[-1.1],[-2.1],[-3.1],[-4.1],[-5.1],[-6.1]])
                    #force_int_m = np.array([[-1.1,-11],[-2.1,-21],[-3.1,-31],[-4.1,-41],[-5.1,-51],[-6.1,-61]])
                    force_int_diff = (force_int_p - force_int_m)
                    #force_int_diff = force_int_p[i,:] - force_int_m[i,:]
                    for s in range(0, deck.dim):
                        if r==s:
                            jacobian[i*deck.dim+r,j*deck.dim+s] = force_int_diff[r] / (2.*eps)
        return jacobian
                            
    def newton_step(self, ysolver, deck, t_n, perturbation_factor, residual):
        jacobian = self.compute_jacobian(ysolver, deck, t_n, perturbation_factor)
        residual = np.reshape(residual,(deck.dim*deck.num_nodes,1))
        
        removeId = []
        for con in deck.conditions:
            if con.type == "Displacement":
                for i in con.id:
                    removeId.append(int((i*deck.dim) + con.direction-1))
        removeId.sort()          

        jacobian = np.delete(jacobian,removeId,0)
        jacobian = np.delete(jacobian,removeId,1)
        residual = np.delete(residual,removeId,0)
        
        delta_y = linalg.solve(jacobian, -residual)
       
        mask = np.ones((deck.num_nodes * deck.dim), dtype=bool)
        mask[removeId] = False
        
        result = np.zeros((deck.num_nodes * deck.dim),dtype=np.float64)
        i = 0
        j = 0
        for m in mask:
            if m == True:
                result[int(i)] = delta_y[int(j)]
                j+= 1
            i += 1
        return np.reshape(result, (deck.num_nodes,deck.dim))

    # This function solves the problem at each time step, using the previous
    # time step solution as an initial guess.
    # This function calls the compute_residual function
    def quasi_static_solver(self, ysolver, deck):
        for t_n in range(1, deck.time_steps):
            res = float('inf')
            iteration = 1
            residual = self.compute_residual(ysolver, deck, t_n)
            res = linalg.norm(residual)
            while res >= deck.solver_tolerance and iteration <= deck.solver_max_it :
                residual = self.compute_residual(ysolver, deck, t_n)
                res = linalg.norm(residual)
                if iteration == deck.solver_max_it:
                    print "Warning: Solver reached limit of " + str(deck.solver_max_it) + " iterations"
                    #sys.exit(1)  
                if res > deck.solver_tolerance:
                    delta_y = self.newton_step(ysolver,deck, t_n, deck.solver_perturbation, residual)
                    ysolver += delta_y
                    residual = self.compute_residual(ysolver, deck, t_n)
                    res = linalg.norm(residual)
                    iteration += 1  
            self.y[:,:,t_n] = ysolver
            print "t_n:" , t_n , "res:" , res , "Iteration #",iteration-1

    # Records the force vector at each time step
    def update_force_data(self, mat_class, t_n):
        self.force_int[:,:, t_n] = mat_class.f_int

    # Records the ext_state vector at each time step
    def update_ext_state_data(self, mat_class, t_n):
        self.ext[:, :, t_n] = mat_class.e

    # Records the ext_state_visco vector at each time step
    def update_ext_state_visco_data(self, mat_class, t_n):
        self.ext_visco[:, :, :, t_n] = mat_class.e_visco

    def strain_calculation(self, id_Node_1, id_Node_2, deck):
        strain = np.zeros( ( deck.time_steps ),dtype=np.float64 )
        for t_n in range(1, deck.time_steps):
            actual = np.linalg.norm(self.y[id_Node_2,:,t_n] - self.y[id_Node_1,:,t_n])
            initial = np.linalg.norm(deck.geometry.nodes[id_Node_2,:] - deck.geometry.nodes[id_Node_1,:])
            strain[t_n] = (actual - initial) / initial
        return strain
            
