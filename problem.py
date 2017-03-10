# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 12:50:28 2015

@author: ilyass.tabiai@polymtl.ca
@author: rolland.delorme@polymtl.ca
@author: patrick.diehl@polymtl.ca
"""

import logging
import numpy as np
import scipy.optimize
import random
#import csv
#import os
#import sys
#import math

logger = logging.getLogger(__name__)

class PD_problem():
    
    def __init__(self, PD_deck):
        # Import initial data
        self.len_x = PD_deck.num_nodes_x
        self.b = np.zeros( ( self.len_x, PD_deck.time_steps) )
        self.compute_b(PD_deck)
        self.compute_horizon(PD_deck)
        self.generate_neighborhood_matrix(PD_deck, PD_deck.geometry.pos_x)
        self.y = np.zeros( ( self.len_x, PD_deck.time_steps) )
        self.y[:,0] = PD_deck.geometry.pos_x
        #self.u = np.zeros( (self.len_x, PD_deck.time_steps ) )
        self.strain = np.zeros( ( PD_deck.time_steps ) )
        self.forces = np.zeros( ( self.len_x, PD_deck.time_steps ) )
        self.ext = np.zeros( ( self.len_x, self.len_x, PD_deck.time_steps ) )
        
        if PD_deck.material_type == "Elastic":
            self.Modulus = PD_deck.e_modulus
        elif PD_deck.material_type == "Viscoelastic":
            self.Relax_Modulus = PD_deck.relax_modulus
            self.Relax_Time = PD_deck.relax_time
            self.ext_visco = np.zeros( ( self.len_x, self.len_x, len(self.Relax_Time), PD_deck.time_steps) ) 
        else:
            logger.error("Error in problem.py: Material type unknown, please use Elastic or Viscoelastic.")      
        
        #self.experimental_nodes = 3
        #self.exp_displacement = np.zeros((PD_deck.time_steps - 1, self.experimental_nodes))
        #self.exp_times = np.zeros((PD_deck.time_steps - 1))
        #self.exp_init_positions = np.zeros(self.experimental_nodes)
        #self.energy = np.zeros( (self.len_x, PD_deck.time_steps ) )

    # Creates a loading vector b which describes the force or displacement applied on each node
    # at any time step
    def compute_b(self, PD_deck):       
        #Build  matrix b[row = node, column = time]
        b = np.zeros( ( self.len_x, PD_deck.time_steps) )
        for t_n in range(1, PD_deck.time_steps): 
            for con in PD_deck.conditions:
                if con.type == "Force":
                    #Madenci approach
                    for x_i in con.id:
                        b[int(x_i), t_n] = self.ramp_loading( PD_deck, t_n , con )
            self.b = b
            #print "b =" , self.b
        
    # Provide the loading shape to use to compute the loading vector b
    def ramp_loading(self, PD_deck, t_n, con):     
        Time_t = PD_deck.delta_t*(t_n)
        if PD_deck.shape_type == "Ramp":
            if con.type == "Force":                                      
                if Time_t <= PD_deck.shape_values[0]:
                    result = (con.force_density*Time_t)/PD_deck.shape_values[0]
                    return result
                elif Time_t > PD_deck.shape_values[0] and Time_t <= PD_deck.shape_values[1]:
                    result = con.force_density
                    return result
                elif Time_t > PD_deck.shape_values[1] and Time_t <= PD_deck.shape_values[2]: 
                    result = con.force_density - con.force_density*(Time_t - PD_deck.shape_values[1])/(PD_deck.shape_values[2] - PD_deck.shape_values[1])
                    return result
                else:
                    result = 0
                    return result
            elif con.type == "Displacement":
                logger.error("Ramp loading in displacement: coming soon. See compute_u_load(PD_deck)")
            else:
                logger.error("Error in problem.py: Type of BC unknown, please use Force or Displacement.")
        else:
            logger.error("Error in problem.py: Shape of BC unknown, please use Ramp.")
                   
    # Computes the horizon        
    def compute_horizon(self, PD_deck):
        #Be sure that points are IN the horizon
        safety_factor = 1.001
        self.Horizon = PD_deck.horizon_factor_m_value*PD_deck.delta_x*safety_factor
        #print "Horizon =" , self.Horizon

    # Returns a list of addresses of the neighbors of a point x_i
    def get_index_x_family(self, x_i):
        return (np.where(self.family[x_i] == 1))[0]

    # Generates matrix neighborhood
    def generate_neighborhood_matrix(self, PD_deck, x):
        self.family = np.zeros((self.len_x, self.len_x))
        for x_i in range(0, len(x)):
            for x_p in range(0, len(x)):
                if x_p == x_i:
                    pass
                elif np.absolute(x_i - x_p) <= PD_deck.horizon_factor_m_value:
                    self.family[x_i][x_p] = 1
                else:
                    pass
        #print self.family
                                                              
    # Computes the shape tensor (here a scalar) for each node
    def compute_m(self, y):
        M = np.zeros((self.len_x, self.len_x))
        for x_i in range(0, self.len_x):
            index_x_family = self.get_index_x_family(x_i)
            for x_p in index_x_family:
                M[x_i, x_p] = (y[x_p] - y[x_i]) / np.absolute(y[x_p] - y[x_i])
        return M

    # Computes the weights for each PD node
    def weighted_function(self, PD_deck, x, x_i):
        index_x_family = self.get_index_x_family(x_i)
        result = 0
        for x_p in index_x_family:
            result = result + PD_deck.influence_function * (x[x_p] - x[x_i])**2 * PD_deck.geometry.volumes[x_p]
        return result

    # Computes the residual vector used in the quasi_static_solver function
    def compute_residual(self, y, PD_deck, t_n):
        residual = np.zeros( ( self.len_x ) )
        for con in PD_deck.conditions:
            if con.type == "Displacement": 
                for id in con.id:
                    y[int(id)] = PD_deck.geometry.pos_x[int(id)] + con.value
        # Choice of the material class
        if PD_deck.material_type == "Elastic":
            from materials.elastic import elastic_material
            variables = elastic_material( PD_deck, self, y )
            self.update_force_data(variables, t_n)
            self.update_ext_state_data(variables, t_n)
        elif PD_deck.material_type == "Viscoelastic":
            from materials.viscoelastic import viscoelastic_material
            variables = viscoelastic_material( PD_deck, self, y, t_n)
            self.update_force_data(variables, t_n)
            self.update_ext_state_data(variables, t_n)
            self.update_ext_state_visco_data(variables, t_n)
        else:
            logger.error("Error in problem.py: Material type unknown, please use Elastic or Viscoelastic.")            
        for x_i in range(0,self.len_x):
            found = False
            for con in PD_deck.conditions:
                if con.type == "Displacement":  
                    if x_i in con.id:
                        found = True
            if found == False:
                residual[x_i] = variables.Ts[x_i] + self.b[x_i, t_n]
        #print residual
        return residual

    # This functtion solves the problem at each time step, using the previous
    # time step solution as an initial guess.
    # This function calls the compute_residual function
    def quasi_static_solver(self, y, PD_deck):
        
        for t_n in range(1, PD_deck.time_steps):
            solver = scipy.optimize.root(self.compute_residual, y, args=(PD_deck, t_n), method='krylov',jac=None,tol=1.0e-12,callback=None,options={'maxiter':1000,'xtol':1.0e-12,'xatol':1.0e-12,'ftol':1.0e-12})
            self.y[:, t_n] = solver.x
            y = self.random_initial_guess(solver.x, PD_deck)
            if solver.success == "False":
                logger.warning("Convergence could not be reached.")
            else:
                logger.info( t_n, solver.success )
            #print y
        return solver

    # Records the force vector at each time step
    def update_force_data(self, variables, t_n):    
        self.forces[:, t_n] = variables.Ts
        
    # Records the ext_state vector at each time step
    def update_ext_state_data(self, variables, t_n):    
        self.ext[:, :, t_n] = variables.e  
    
    # Records the ext_state_visco vector at each time step
    def update_ext_state_visco_data(self, variables, t_n):    
        self.ext_visco[:, :, :, t_n] = variables.e_visco  

    # Initial guess
    def random_initial_guess(self, z, PD_deck):
        #Do not forget to do this for each direction, not only x
        y = np.zeros((self.len_x))
        y = z + 0.25 * random.uniform(-1, 1) * PD_deck.delta_x
        return y

    def strain_center_bar(self, PD_deck):
        Mid_Node_1 = int(PD_deck.num_nodes_x/2)-1
        #print "Mid_Node_1 =" , Mid_Node_1
        Mid_Node_2 = int(PD_deck.num_nodes_x/2)+1
        #print "Mid_Node_2 =" , Mid_Node_2
        for t_n in range(1, PD_deck.time_steps):
            self.strain[t_n] = (np.absolute(self.y[Mid_Node_2,t_n] - self.y[Mid_Node_1,t_n]) - np.absolute(PD_deck.geometry.pos_x[Mid_Node_2] - PD_deck.geometry.pos_x[Mid_Node_1])) / np.absolute(PD_deck.geometry.pos_x[Mid_Node_2] - PD_deck.geometry.pos_x[Mid_Node_1])


#    # Records the force vector at each time step
#    def update_energy_data(self, variables, t_n):
#        self.energy = self.strain_energy_from_force
#
#    # Records the ext_state vector at each time step
#    def update_displacements(self, t_n):
#        self.u[:, t_n] = self.y[:, t_n] - self.y[:, 0]
#
#    # Computes the strain energy density from Ts x u
#    def strain_energy_from_force(self, PD_deck):
#        energy = np.zeros( (PD_deck.time_steps, self.len_x) )
#        for t_n in range(0, PD_deck.time_steps ):   
#            for x_i in range(0, PD_deck.num_nodes_x):               
#                energy[t_n , x_i] = abs(self.forces[x_i, t_n]) * abs(self.u[x_i, t_n]) * PD_deck.Volume
#                #abs(self.forces[x_i, t_n]), abs(self.u[x_i, t_n]), energy[t_n , x_i]
#        #print "ENERGY:",
#        #print energy
#        self.strain_energy_from_force = energy         
#    
#    # Computes the strain energy using the formula in the PMB
#    def strain_energy_bond_based(self, PD_deck):
#        energy = np.zeros((self.len_x, PD_deck.time_steps))
#        for t_n in range(0, PD_deck.time_steps):
#            for x_i in range(0, PD_deck.num_nodes_x):
#                index_x_family = self.get_index_x_family(x_i)
#                modulus = PD_deck.get_elastic_material_properties()
#                for x_p in index_x_family:
#                    #Silling-Askari2005, Eq17
#                    stretch = (abs(self.u[x_p, t_n] - self.u[x_i, t_n]) - abs(
#                        PD_deck.geometry.pos_x[x_p] - PD_deck.geometry.pos_x[x_i])) / abs(PD_deck.geometry.pos_x[x_p] - PD_deck.geometry.pos_x[x_i])
#                    #Silling-Askari2005, Eq22
#                    energy[x_i, t_n] = (
#                        math.pi * modulus * math.pow(stretch, 2) * math.pow(self.Horizon, 4)) / 4.
#        self.strain_energy = energy
