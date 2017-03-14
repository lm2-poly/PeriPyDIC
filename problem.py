#-*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrick.diehl@polymtl.ca

import logging
import numpy as np
import scipy.optimize
import random
import util.neighbor

logger = logging.getLogger(__name__)

class PD_problem():
    
    def __init__(self, deck):
        # Import initial data
        self.len_x = deck.num_nodes_x
        self.b = np.zeros( ( self.len_x, deck.time_steps) )
        self.compute_b(deck)
        self.neighbors = util.neighbor.NeighborSearch(deck)
        self.y = np.zeros( ( self.len_x, deck.time_steps) )
        self.y[:,0] = deck.geometry.pos_x
        #self.u = np.zeros( (self.len_x, deck.time_steps ) )
        self.strain = np.zeros( ( deck.time_steps ) )
        self.forces = np.zeros( ( self.len_x, deck.time_steps ) )
        self.ext = np.zeros( ( self.len_x, self.len_x, deck.time_steps ) )
        
        if deck.material_type == "Elastic":
            self.Modulus = deck.e_modulus
        elif deck.material_type == "Viscoelastic":
            self.Relax_Modulus = deck.relax_modulus
            self.Relax_Time = deck.relax_time
            self.ext_visco = np.zeros( ( self.len_x, self.len_x, len(self.Relax_Time), deck.time_steps) ) 
        else:
            logger.error("Error in problem.py: Material type unknown, please use Elastic or Viscoelastic.")      
        
        #self.experimental_nodes = 3
        #self.exp_displacement = np.zeros((deck.time_steps - 1, self.experimental_nodes))
        #self.exp_times = np.zeros((deck.time_steps - 1))
        #self.exp_init_positions = np.zeros(self.experimental_nodes)
        #self.energy = np.zeros( (self.len_x, deck.time_steps ) )

    # Creates a loading vector b which describes the force or displacement applied on each node
    # at any time step
    def compute_b(self, deck):       
        #Build  matrix b[row = node, column = time]
        b = np.zeros( ( self.len_x, deck.time_steps) )
        for t_n in range(1, deck.time_steps): 
            for con in deck.conditions:
                if con.type == "Force":
                    #Madenci approach
                    for x_i in con.id:
                        b[int(x_i), t_n] = self.ramp_loading( deck, t_n , con )
            self.b = b
            #print "b =" , self.b
        
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
                   
    # Computes the shape tensor (here a scalar) for each node
    def compute_m(self, y):
        M = np.zeros((self.len_x, self.len_x))
        for x_i in range(0, self.len_x):
            index_x_family = self.neighbors.get_index_x_family(x_i)
            for x_p in index_x_family:
                M[x_i, x_p] = (y[x_p] - y[x_i]) / np.absolute(y[x_p] - y[x_i])
        return M

    # Computes the weights for each PD node
    def weighted_function(self, deck, x, x_i):
        index_x_family = self.neighbors.get_index_x_family(x_i)
        result = 0
        for x_p in index_x_family:
            result = result + deck.influence_function * (x[x_p] - x[x_i])**2 * deck.geometry.volumes[x_p]
        return result

    # Computes the residual vector used in the quasi_static_solver function
    def compute_residual(self, y, deck, t_n):
        residual = np.zeros( ( self.len_x ) )
        for con in deck.conditions:
            if con.type == "Displacement": 
                for id in con.id:
                    y[int(id)] = deck.geometry.pos_x[int(id)] + con.value
        # Choice of the material class
        if deck.material_type == "Elastic":
            from materials.elastic import Elastic_material
            variables = Elastic_material( deck, self, y )
            self.update_force_data(variables, t_n)
            self.update_ext_state_data(variables, t_n)
        elif deck.material_type == "Viscoelastic":
            from materials.viscoelastic import Viscoelastic_material
            variables = Viscoelastic_material( deck, self, y, t_n)
            self.update_force_data(variables, t_n)
            self.update_ext_state_data(variables, t_n)
            self.update_ext_state_visco_data(variables, t_n)
        else:
            logger.error("Error in problem.py: Material type unknown, please use Elastic or Viscoelastic.")            
        for x_i in range(0,self.len_x):
            found = False
            for con in deck.conditions:
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
    def quasi_static_solver(self, y, deck):
        
        for t_n in range(1, deck.time_steps):
            solver = scipy.optimize.root(self.compute_residual, y, args=(deck, t_n), method=deck.solver_type,jac=None,tol=deck.solver_tolerance,callback=None,options={'maxiter':1000,'xtol':1.0e-12,'xatol':1.0e-12,'ftol':1.0e-12})
            self.y[:, t_n] = solver.x
            y = self.random_initial_guess(solver.x, deck)
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
    def random_initial_guess(self, z, deck):
        #Do not forget to do this for each direction, not only x
        y = np.zeros((self.len_x))
        y = z + 0.25 * random.uniform(-1, 1) * deck.delta_x
        return y

    def strain_center_bar(self, deck):
        Mid_Node_1 = int(deck.num_nodes_x/2)-1
        #print "Mid_Node_1 =" , Mid_Node_1
        Mid_Node_2 = int(deck.num_nodes_x/2)+1
        #print "Mid_Node_2 =" , Mid_Node_2
        for t_n in range(1, deck.time_steps):
            self.strain[t_n] = (np.absolute(self.y[Mid_Node_2,t_n] - self.y[Mid_Node_1,t_n]) - np.absolute(deck.geometry.pos_x[Mid_Node_2] - deck.geometry.pos_x[Mid_Node_1])) / np.absolute(deck.geometry.pos_x[Mid_Node_2] - deck.geometry.pos_x[Mid_Node_1])


#    # Records the force vector at each time step
#    def update_energy_data(self, variables, t_n):
#        self.energy = self.strain_energy_from_force
#
#    # Records the ext_state vector at each time step
#    def update_displacements(self, t_n):
#        self.u[:, t_n] = self.y[:, t_n] - self.y[:, 0]
#
#    # Computes the strain energy density from Ts x u
#    def strain_energy_from_force(self, deck):
#        energy = np.zeros( (deck.time_steps, self.len_x) )
#        for t_n in range(0, deck.time_steps ):   
#            for x_i in range(0, deck.num_nodes_x):               
#                energy[t_n , x_i] = abs(self.forces[x_i, t_n]) * abs(self.u[x_i, t_n]) * deck.Volume
#                #abs(self.forces[x_i, t_n]), abs(self.u[x_i, t_n]), energy[t_n , x_i]
#        #print "ENERGY:",
#        #print energy
#        self.strain_energy_from_force = energy         
#    
#    # Computes the strain energy using the formula in the PMB
#    def strain_energy_bond_based(self, deck):
#        energy = np.zeros((self.len_x, deck.time_steps))
#        for t_n in range(0, deck.time_steps):
#            for x_i in range(0, deck.num_nodes_x):
#                index_x_family = self.get_index_x_family(x_i)
#                modulus = deck.get_elastic_material_properties()
#                for x_p in index_x_family:
#                    #Silling-Askari2005, Eq17
#                    stretch = (abs(self.u[x_p, t_n] - self.u[x_i, t_n]) - abs(
#                        deck.geometry.pos_x[x_p] - deck.geometry.pos_x[x_i])) / abs(deck.geometry.pos_x[x_p] - deck.geometry.pos_x[x_i])
#                    #Silling-Askari2005, Eq22
#                    energy[x_i, t_n] = (
#                        math.pi * modulus * math.pow(stretch, 2) * math.pow(self.Horizon, 4)) / 4.
#        self.strain_energy = energy
