#-*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrick.diehl@polymtl.ca

import logging
import numpy as np
from numpy import linalg
import scipy.optimize
import random
import util.neighbor
np.set_printoptions(threshold='nan')
import sys

logger = logging.getLogger(__name__)

class PD_problem():

    def __init__(self, deck):
        # Import initial data

        # NeighborSearch
        self.neighbors = util.neighbor.NeighborSearch(deck)

        #self.b = np.zeros( ( deck.dim * deck.num_nodes, deck.time_steps) )
        #Compute the force density as a 3D array with the last dimension for the time steps
        self.compute_b(deck)

        #Compute the weighted volume for each node in a vector.
        self.weighted_function(deck)

        #self.y = np.zeros( ( deck.dim * deck.num_nodes, deck.time_steps) )
        #Actual position as a 3D array with the last dimension for the time steps
        self.y = np.zeros((deck.num_nodes, deck.dim, deck.time_steps))
        self.y[:,:,0] = deck.geometry.nodes[:,:]

        self.strain = np.zeros( ( deck.time_steps ) )
        self.forces = np.zeros((deck.num_nodes, deck.dim, deck.time_steps))
        self.ext = np.zeros( ( deck.num_nodes, deck.num_nodes, deck.time_steps ) )

    # Creates a loading vector b which describes the force or displacement applied on each node
    # at any time step
    def compute_b(self, deck):
        #Force density as a 3D array with the last dimension for the time steps
        self.b = np.zeros((deck.num_nodes, deck.dim, deck.time_steps))
        for t_n in range(1, deck.time_steps):
            for con in deck.conditions:
                if con.type == "Force":
                    #Madenci approach
                    for i in con.id:
                        # x direction
                        if con.direction == 1:
                            self.b[int(i), 0, t_n] = self.ramp_loading( deck, t_n , con )
                        # y direction
                        if con.direction == 2:
                            self.b[int(i),  1 , t_n] = self.ramp_loading( deck, t_n , con )
                        # z direction
                        if con.direction == 3:
                            self.b[int(i), 2, t_n] = self.ramp_loading( deck, t_n , con )
        #print self.b

    # Provide the loading shape to use to compute the loading vector b
    def ramp_loading(self, deck, t_n, con):
        Time_t = deck.delta_t*(t_n)
        if deck.shape_type == "Ramp":
            if con.type == "Force":
                if Time_t <= deck.shape_values[0]:
                    result = (con.force_density*Time_t)/deck.shape_values[0]
                    return result
                elif Time_t > deck.shape_values[0] and Time_t <= deck.shape_values[1]:
                    result = con.force_density
                    return result
                elif Time_t > deck.shape_values[1] and Time_t <= deck.shape_values[2]:
                    result = con.force_density - con.force_density*(Time_t - deck.shape_values[1])/(deck.shape_values[2] - deck.shape_values[1])
                    return result
                else:
                    result = 0
                    return result
            elif con.type == "Displacement":
                logger.error("Ramp loading in displacement: coming soon. See compute_u_load(deck)")
            else:
                logger.error("Error in problem.py: Type of BC unknown, please use Force or Displacement.")
        else:
            logger.error("Error in problem.py: Shape of BC unknown, please use Ramp.")

    # Computes the weights for each PD node
    def weighted_function(self, deck):
        self.weighted_volume = np.zeros((deck.num_nodes))
        for i in range(0, deck.num_nodes):
            index_x_family = self.neighbors.get_index_x_family(i)
            for p in index_x_family:
                X = deck.geometry.nodes[p,:] - deck.geometry.nodes[i,:]
                self.weighted_volume[i] += deck.influence_function * (np.linalg.norm(X))**2 * deck.geometry.volumes[p]

    # Computes the residual vector used in the quasi_static_solver function
    def compute_residual(self, ysolver, deck, t_n):
        residual = np.zeros((deck.num_nodes, deck.dim))
        for con in deck.conditions:
            if con.type == "Displacement":
                for id_node in con.id:
                    # x direction
                    if con.direction == 1:
                        if deck.dim == 1:
                            ysolver[int(id_node),0] = deck.geometry.nodes[int(id_node),0] + con.value
                    # y direction
                    if con.direction == 2:
                        ysolver[int(id_node),1] = deck.geometry.nodes[int(id_node),1] + con.value
                    # z direction
                    if con.direction == 3:
                        ysolver[int(id_node),2] = deck.geometry.nodes[int(id_node),2] + con.value

        # Choice of the material class
        if deck.material_type == "Elastic":
            from materials.elastic import Elastic_material
            variables = Elastic_material( deck, self, ysolver )
            self.update_force_data(variables, t_n)
            self.update_ext_state_data(variables, t_n)
        elif deck.material_type == "Viscoelastic":
            from materials.viscoelastic import Viscoelastic_material
            variables = Viscoelastic_material( deck, self, ysolver, t_n)
            self.update_force_data(variables, t_n)
            self.update_ext_state_data(variables, t_n)
            self.update_ext_state_visco_data(variables, t_n)
        else:
            logger.error("Error in problem.py: Material type unknown, please use Elastic or Viscoelastic.")
        for i in range(0,deck.num_nodes):
            found = False
            for con in deck.conditions:
                if con.type == "Displacement":
                    if i in con.id:
                        found = True
            if found == False:
                residual[i,:] = variables.F[i,:] + self.b[i,:, t_n]
        #print "dila", variables.pd_dilatation(deck, self, deck.geometry.nodes, variables.e, 16)
        #print residual
        return residual

    def compute_jacobian(self, deck, y_actual, var_actual ):
        y_perturb = y_actual.copy()
        epsilon = 1.0e-6 * deck.delta_X
        y_perturb += epsilon
        
        if deck.material_type == "Elastic":
            from materials.elastic import Elastic_material
            var_perturb = Elastic_material( deck, self, y_perturb )
            
        self.jocabian = np.zeros((deck.num_nodes * deck.dim , deck.num_nodes * deck.dim))
        for m in range(0, deck.num_nodes)
            for k in self.neighbors.get_index_x_family(m):
                for i in range(0, deck.dim):
                    for j in range(0, deck.dim):
                        self.jacobian[m*deck.dim+i,k*deck.dim+j] = (var_pertub.F[k,m] - var_actual.F[k,m]) / epsilon









    # This function solves the problem at each time step, using the previous
    # time step solution as an initial guess.
    # This function calls the compute_residual function
    def quasi_static_solver(self, ysolver, deck):

        for t_n in range(1, deck.time_steps):
            solver = scipy.optimize.root(self.compute_residual, ysolver, args=(deck, t_n), method=deck.solver_type,jac=None,tol=None,callback=None,options={'maxiter':1000, 'ftol':deck.solver_tolerance, 'fatol':deck.solver_tolerance})
            self.y[:,:,t_n] = solver.x[:,:]
            ysolver = self.random_initial_guess(solver.x, deck)
            if solver.success == "False":
                print "Convergence could not be reached."
            else:
                print "Time Step: ", t_n, "Convergence: ", solver.success
            #print ysolver
            #print t_n
        return solver

    # Records the force vector at each time step
    def update_force_data(self, variables, t_n):
        self.forces[:,:, t_n] = variables.F

    # Records the ext_state vector at each time step
    def update_ext_state_data(self, variables, t_n):
        self.ext[:, :, t_n] = variables.e

    # Records the ext_state_visco vector at each time step
    def update_ext_state_visco_data(self, variables, t_n):
        self.ext_visco[:, :, :, t_n] = variables.e_visco

    # Initial guess
    def random_initial_guess(self, z, deck):
        y = np.zeros((deck.num_nodes, deck.dim))
        for i in range(0,deck.num_nodes):
            for j in range(0,deck.dim):
                y[i,j] = z[i,j] + 0.25 * random.uniform(-1, 1) * deck.delta_X
        return y

    def strain_calculation(self, id_Node_1, id_Node_2, deck):
        for t_n in range(1, deck.time_steps):
            actual = np.linalg.norm(self.y[id_Node_2,:,t_n] - self.y[id_Node_1,:,t_n])
            initial = np.linalg.norm(deck.geometry.nodes[id_Node_2,:] - deck.geometry.nodes[id_Node_1,:])
            self.strain[t_n] = (actual - initial) / initial
            
            

