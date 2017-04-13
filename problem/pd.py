#-*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrick.diehl@polymtl.ca

import logging
import numpy as np
import util.neighbor
from scipy import linalg
np.set_printoptions(precision=8, threshold='nan', suppress=True)
#import sys

logger = logging.getLogger(__name__)

## Class to define the peridynamic problem, i.e. applying boundaries conditions to the geomety and solve the problem
class PD_problem():
    
    ## Constructor
    # @param deck The input deck   
    def __init__(self, deck):
       
        ## Family of each node
        self.neighbors = util.neighbor.NeighborSearch(deck)
        
        ## Nodes' positions stored for each time step
        self.y = np.zeros((deck.num_nodes, deck.dim, deck.time_steps),dtype=np.float64)
        self.y[:,:,0] = deck.geometry.nodes[:,:]
        
        ## Global internal force density array storing the force density attached to each node for each time step
        self.force_int = np.zeros((deck.num_nodes, deck.dim, deck.time_steps),dtype=np.float64)
        
        ## Extension array storing the extension at each node between the node and its family                7
        self.ext = np.zeros( ( deck.num_nodes, deck.num_nodes, deck.time_steps ),dtype=np.float64 )

        ## Compute the external force density "b" applied on each node
        self.compute_b(deck)

        # Compute the weighted volume for each node
        self.compute_weighted_volume(deck)

    ## Compute the external force density "b" applied on each node
    # @param deck The input deck
    def compute_b(self, deck):
        ## External force density "b" applied on each node
        self.b = np.zeros((deck.num_nodes, deck.dim, deck.time_steps),dtype=np.float64)
        for t_n in range(1, deck.time_steps):
            for con in deck.conditions:
                if con.type == "Force":
                    if con.shape == "Ramp":
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

    ## Provide the loading shape
    # @param deck The input deck
    # @param t_n Id of the time step
    # @param con Type of loading, "Force" or "Displacement"
    # @param i Id of Node "i"
    def shape_loading(self, deck, t_n, con, i):
        Time_t = deck.delta_t*(t_n)
        if deck.shape_type == "Ramp":
            
            if con.type == "Force":
                value = con.value / deck.geometry.volumes[int(i)]
            if con.type == "Discplacement":
                value = con.value 
                
            if Time_t <= deck.shape_values[0]:
                result = (value*Time_t)/deck.shape_values[0]
                return result
            elif Time_t > deck.shape_values[0] and Time_t <= deck.shape_values[1]:
                result = value
                return result
            elif Time_t > deck.shape_values[1] and Time_t <= deck.shape_values[2]:
                result = value - value*(Time_t - deck.shape_values[1])/(deck.shape_values[2] - deck.shape_values[1])
                return result
            else:
                result = 0
                return result
        else:
            logger.error("Error in problem.py: Shape of BC unknown, please use Ramp.")

    ## Compute the weighted volume for each node
    # @param deck The input deck
    def compute_weighted_volume(self, deck):
        ## Weighted volume for each node
        self.weighted_volume = np.zeros((deck.num_nodes),dtype=np.float64)
        for i in range(0, deck.num_nodes):
            index_x_family = self.neighbors.get_index_x_family(i)
            for p in index_x_family:
                X = deck.geometry.nodes[p,:] - deck.geometry.nodes[i,:]
                self.weighted_volume[i] += deck.influence_function * (np.linalg.norm(X))**2 * deck.geometry.volumes[p]

    ## Provide the internal force density for each node for a given time step t_n
    # @param deck The input deck    
    # @param ysolver Initial guess for Actual nodes' position
    # @param t_n Id of the time step
    # @return Internal force density for each node  
    def internal_force(self, deck, ysolver, t_n):
        # Choice of the material class
        if deck.material_type == "Elastic":
            from materials.elastic import Elastic_material
            ## Data from the material class            
            self.mat_class = Elastic_material( deck, self, ysolver )
            # Data from the material class            
            self.update_force_data(self.mat_class, t_n)
            # Data from the material class            
            self.update_ext_state_data(self.mat_class, t_n)
            
        elif deck.material_type == "Viscoelastic":
            from materials.viscoelastic import Viscoelastic_material
            self.mat_class = Viscoelastic_material( deck, self, ysolver, t_n)
            self.update_force_data(self.mat_class, t_n)
            self.update_ext_state_data(self.mat_class, t_n)
            self.update_ext_state_visco_data(self.mat_class, t_n)
        else:
            logger.error("Error in problem.py: Material type unknown, please use Elastic or Viscoelastic.")
        
        force = self.mat_class.f_int
        #internal_force = np.reshape(internal_force, (deck.num_nodes * deck.dim,-1) ) 
        return force
    
    ## Provide the residual for each node after a solving step for a given time step t_n
    # @param deck The input deck    
    # @param ysolver Initial guess for Actual nodes' position
    # @param t_n Id of the time step
    # @return Residual for each node
    def residual_vector(self, deck, ysolver, t_n):
        residual = np.zeros((deck.num_nodes, deck.dim),dtype=np.float64)
        internal_force = self.internal_force(deck, ysolver, t_n)
        for con in deck.conditions:
            if con.type == "Displacement" and con.shape == "Fixed":
                for id_node in con.id:
                    # x direction
                    if con.direction == 1:
                        ysolver[int(id_node),0] = deck.geometry.nodes[int(id_node),0] + con.value
                    # y direction
                    if con.direction == 2:
                        ysolver[int(id_node),1] = deck.geometry.nodes[int(id_node),1] + con.value
                    # z direction
                    if con.direction == 3:
                        ysolver[int(id_node),2] = deck.geometry.nodes[int(id_node),2] + con.value
                        
            if con.type == "Displacement" and con.shape == "Ramp":
                for i in con.id:
                    # x direction
                    if con.direction == 1:
                        ysolver[int(id_node),0] = self.shape_loading( deck, t_n , con , i )
                    # y direction
                    if con.direction == 2:
                        ysolver[int(id_node),1] = self.shape_loading( deck, t_n , con , i )
                    # z direction
                    if con.direction == 3:
                        ysolver[int(id_node),2] = self.shape_loading( deck, t_n , con , i )   
                         
        for i in range(0,deck.num_nodes):
            found = False
            for con in deck.conditions:
                if con.type == "Displacement":
                    if i in con.id:
                        found = True
            if found == False:  
                residual[i,:] = internal_force[i,:] + self.b[i,:, t_n]               
        return residual

    ## Provide the Jacobian (stiffness) matrix for a given time step t_n for the Newton's method
    # @param deck The input deck    
    # @param ysolver Initial guess for Actual nodes' position
    # @param t_n Id of the time step
    # @param perturbation_factor Magnitude of the perturbation factor
    # @return Jacobian matrix
    def jacobian_matrix(self, deck, ysolver, t_n, perturbation_factor):
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
                    force_int_p = self.internal_force(deck, ysolver + eps_vector, t_n)[i,:]
                    force_int_m = self.internal_force(deck, ysolver - eps_vector, t_n)[i,:]
                    force_int_diff = (force_int_p - force_int_m)
                    for s in range(0, deck.dim):
                        if r==s:
                            jacobian[i*deck.dim+r,j*deck.dim+s] = force_int_diff[r] / (2.*eps)
        #print "Jacobian Matrix Density =", np.count_nonzero(jacobian) / (float(deck.num_nodes) * float(deck.dim))**2 * 100., "%"
        return jacobian

    ## Provide the displacement increment resulting for the Newton's method, for each node for a given time step t_n
    # @param deck The input deck    
    # @param ysolver Initial guess for Actual nodes' position
    # @param t_n Id of the time step
    # @param perturbation_factor Magnitude of the perturbation factor
    # @param residual Residual for each node resulting from a solving step
    # @return Displacement increment for each node                            
    def newton_step(self, deck, ysolver, t_n, perturbation_factor, residual):
        jacobian = self.jacobian_matrix(deck, ysolver, t_n, perturbation_factor)
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
        
        #delta_y = linalg.solve(jacobian, -residual, check_finite = "False", assume_a = "sym" )
        delta_y = linalg.solve(jacobian, -residual, check_finite = "False")
       
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

    ## Solve the peridynamic problem at each time step using the Newton's method to obtain the actual nodes' position
    # @param deck The input deck    
    # @param ysolver Initial guess for Actual nodes' position
    def quasi_static_solver(self, deck, ysolver):
        for t_n in range(1, deck.time_steps):
            res = float('inf')
            iteration = 1
            residual = self.residual_vector(deck, ysolver, t_n)
            res = linalg.norm(residual)
            while res >= deck.solver_tolerance and iteration <= deck.solver_max_it :             
                print "iteration", iteration                
                if iteration == deck.solver_max_it:
                    print "Warning: Solver reached limit of " + str(deck.solver_max_it) + " iterations"  
                delta_y = self.newton_step(deck, ysolver, t_n, deck.solver_perturbation, residual)
                ysolver += delta_y
                residual = self.residual_vector(deck, ysolver, t_n)
                res = linalg.norm(residual)
                iteration += 1  
            self.y[:,:,t_n] = ysolver
            print "t_n:" , t_n , "res:" , res , "Iteration #",iteration-1

    ## Store the internal force density for each node at each time step
    # @param mat_class Data from the material class
    # @param t_n Id of the time step
    def update_force_data(self, mat_class, t_n):
        ## Global internal force density array storing the force density attached to each node for each time step        
        self.force_int[:,:, t_n] = mat_class.f_int

    ## Store the extension for each node between itself and its family
    # @param mat_class Data from the material class
    # @param t_n Id of the time step 
    def update_ext_state_data(self, mat_class, t_n):
        ## Extension array storing the extension at each node between the node and its family                        
        self.ext[:, :, t_n] = mat_class.e

    ## Store the viscous extension for each node between itself and its family
    # @param mat_class Data from the material class
    # @param t_n Id of the time step
    def update_ext_state_visco_data(self, mat_class, t_n):
        ## Viscous extension array storing the extension at each node between the node and its family        
        self.ext_visco[:, :, :, t_n] = mat_class.e_visco

    ## Provide the strain between 2 nodes
    # @param deck The input deck  
    # @param id_Node_1 Id of the 1st node
    # @param id_Node_2 Id of the 2nd node
    def strain_calculation(self, deck, id_Node_1, id_Node_2):
        strain = np.zeros( ( deck.time_steps ),dtype=np.float64 )
        for t_n in range(1, deck.time_steps):
            actual = np.linalg.norm(self.y[id_Node_2,:,t_n] - self.y[id_Node_1,:,t_n])
            initial = np.linalg.norm(deck.geometry.nodes[id_Node_2,:] - deck.geometry.nodes[id_Node_1,:])
            strain[t_n] = (actual - initial) / initial
        return strain
            
