# -*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrick.diehl@polymtl.ca
import numpy as np
from scipy import linalg
np.set_printoptions(precision=8, threshold='nan')
import sys

## Class to compute the wll-known strain and stress tensors defined in the classical continuum mechanics
class CCM_calcul():

    ## Constructor
    # @param deck The input deck
    # @param data_solver Data from the peridynamic problem/solving class
    def __init__(self, deck, data_solver):
        
        ## Nodes' initial positions
        self.x = deck.geometry.nodes
        
        ## Nodes' positions stored for each time step
        self.y = data_solver.y
        
        ## Global internal force density vector storing the force density attached to each node for each time step
        self.force_int = data_solver.force_int
        
        ## Extension matrix storing the extension at each node between the node and its family        
        self.ext = data_solver.ext
        
        ## Dimension of the data_solver (1D, 2D or 3D)
        self.dim = deck.dim
        
        ## Amount of nodes in the data_solver
        self.num_nodes = deck.num_nodes
        
        ## Amount of time step
        self.time_steps = deck.time_steps
        
        ## Volume related to each node
        self.node_volumes = deck.geometry.volumes        
        
        ## Influence function
        self.influence_function = deck.influence_function
        
        ## Golbal strain tensor storing the strain tensor for each node at each time step      
        self.global_strain_tensor(data_solver)

#        self.compute_u_displacement()
#        self.global_stress_tensor(data_solver)
           
    ## Return the image of (xi - xp) under the reference position vector state X 
    # @param data_solver Data from the peridynamic problem/solving class
    # @param i Id of node #i
    # @param p Id of node #p with node #i family 
    # @return the image 
    def X_vector_state(self, data_solver, i, p):
        image = self.x[p,:] - self.x[i,:]
        image = np.reshape(image,(self.dim,1))
        return image

    ## Return the image of (xi - xp) under the deformation vector state Y
    # @param data_solver Data from the peridynamic problem/solving class
    # @param i Id of Node #i
    # @param p Id of Node #p with Node #i family
    # @param t_n Id of the time step  
    # @return the image  
    def Y_vector_state(self, data_solver, i, p, t_n):
        image = self.y[p,:,t_n] - self.y[i,:,t_n]
        image = np.reshape(image,(self.dim,1))
        return image
        
    ## Return the shape tensor K related to node i
    # @param data_solver Data from the peridynamic problem/solving class
    # @param i Id of Node #i
    def K_shape_tensor(self, data_solver, i):
        K = np.zeros((self.dim, self.dim),dtype=np.float64)
        index_x_family = data_solver.neighbors.get_index_x_family(i)
        for p in index_x_family:
            X = self.X_vector_state(data_solver, i, p)
            K += self.influence_function * np.dot(X,X.T) * self.node_volumes[p]
        return K

    ## Return the deformation gradient tensor epsilon related to Node #i
    # @param data_solver Data from the peridynamic problem/solving class
    # @param i Id of Node #i
    # @param t_n Id of the time step  
    def deformation_gradient(self, data_solver, i, t_n):
        tmp = np.zeros((self.dim, self.dim),dtype=np.float64)       
        index_x_family = data_solver.neighbors.get_index_x_family(i)
        for p in index_x_family:
            Y = self.Y_vector_state(data_solver, i, p, t_n)            
            X = self.X_vector_state(data_solver, i, p)
            tmp += self.influence_function * np.dot(Y,X.T) * self.node_volumes[p]
        deformation = np.dot(tmp, linalg.inv(self.K_shape_tensor(data_solver, i)))
        return deformation
        
    ## Return the strain tensor related to node i
    # @param data_solver Data from the peridynamic problem/solving class
    # @param i Id of Node #i
    # @param t_n Id of the time step 
    def strain_tensor(self, data_solver, i, t_n):
        F = self.deformation_gradient(data_solver, i, t_n)
        strain = (F + F.T)/2 - np.identity(self.dim, dtype=np.float64)
        return strain

    ## Compute the strain tensor for each node
    # @param data_solver Data from the peridynamic problem/solving class
    def global_strain_tensor(self, data_solver):
        ## Golbal strain tensor storing the strain tensor for each node at each time step
        self.global_strain = np.zeros((self.num_nodes*self.dim, self.dim, self.time_steps),dtype=np.float64)        
        for t_n in range(1, self.time_steps):
            for i in range(0, self.num_nodes):
                for j in range(0, self.dim):
                    for r in range(0, self.dim):                    
                        self.global_strain[i*self.dim+r,j,t_n] = self.strain_tensor(data_solver, i, t_n)[r,j]

    ## Compute the displacement for each node at each time step
    def compute_u_displacement(self):
        ## Displacement vector between two consecutives time steps for each node
        self.u = np.zeros((self.num_nodes, self.dim, self.time_steps),dtype=np.float64)
        for t_n in range(1, self.time_steps):
            for i in range(0, self.num_nodes):
                self.u[i,:,t_n] = self.y[i,:,t_n] - self.y[i,:,t_n-1]

#    # Compute the stiffness tensor
#    def compute_C_stiffness_tensor(self, deck):      
#        delta = np.identity(self.dim, dtype=np.float64)
#        
#        J = np.zeros((self.dim, self.dim, self.dim, self.dim),dtype=np.float64)
#        for i in range(0, self.dim):
#            for j in range(0, self.dim):
#                for k in range(0, self.dim):
#                    for l in range(0, self.dim):
#                        J[i,j,k,l] = (1./3.) * delta[i,j]*delta[k,l]
#
#        I = np.zeros((self.dim, self.dim, self.dim, self.dim),dtype=np.float64)
#        for i in range(0, self.dim):
#            for j in range(0, self.dim):
#                for k in range(0, self.dim):
#                    for l in range(0, self.dim):
#                        I[i,j,k,l] = (1./2.) * ( delta[i,k]*delta[j,l] + delta[i,l]*delta[j,k] )
#
#        K = np.zeros((self.dim, self.dim, self.dim, self.dim),dtype=np.float64)
#        for i in range(0, self.dim):
#            for j in range(0, self.dim):
#                for k in range(0, self.dim):
#                    for l in range(0, self.dim):
#                        K[i,j,k,l] = I[i,j,k,l] - J[i,j,k,l]
#        
#        if self.dim ==1:
#            bulk_modulus = deck.young_modulus / 3.
#            shear_modulus = deck.young_modulus /2.
#        
#        if self.dim >= 2:
#            bulk_modulus = deck.bulk_modulus
#            shear_modulus = deck.shear_modulus
#            
#        self.C = 3 * bulk_modulus * J + 2 * shear_modulus * K
#        
#    # Return the stress tensor related to node i
#    def stress_tensor(self, data_solver, i, t_n):
#        stress = np.zeros((self.dim, self.dim),dtype=np.float64)
#        for h in range(0, self.dim):
#            for j in range(0, self.dim):
#                for k in range(0, self.dim):
#                    for l in range(0, self.dim):
#                        stress[h,j] += self.C[h,j,k,l] * self.strain_tensor(data_solver, i, t_n)[k,l]
#                        
#        return stress
#
#    # Return the stress tensor for each node
#    def global_stress_tensor(self, data_solver):
#        self.global_stress = np.zeros((self.num_nodes*self.dim, self.dim, self.time_steps),dtype=np.float64)        
#        for t_n in range(1, self.time_steps):
#            for i in range(0, self.num_nodes):
#                for j in range(0, self.dim):
#                    for r in range(0, self.dim):                    
#                        self.global_stress[i*self.dim+r,j,t_n] = self.stress_tensor(data_solver, i, t_n)[r,j]            