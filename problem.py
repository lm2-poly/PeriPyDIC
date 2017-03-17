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
        self.len_x = deck.dim *  deck.num_nodes
        self.b = np.zeros( ( self.len_x, deck.time_steps) )
        self.compute_b(deck)
        self.neighbors = util.neighbor.NeighborSearch(deck)
        self.y = np.zeros( ( self.len_x, deck.time_steps) )
        if deck.dim == 1:
            self.y[:,0] = deck.geometry.nodes[:,0]
        if deck.dim == 2:
            coordinates = [row[0] for row in deck.geometry.nodes ] + [row[1] for row in deck.geometry.nodes ]
            self.y[:,0] = coordinates
        if deck.dim == 3:
            coordinates = [row[0] for row in deck.geometry.nodes ] + [row[1] for row in deck.geometry.nodes ] + [row[2] for row in deck.geometry.nodes ]
            self.y[:,0] = coordinates
        #self.u = np.zeros( (self.len_x, deck.time_steps ) )
        self.strain = np.zeros( ( deck.time_steps ) )
        self.forces = np.zeros( ( self.len_x, deck.time_steps ) )
        self.ext = np.zeros( ( deck.num_nodes, deck.num_nodes, deck.time_steps ) )

        print deck.geometry.nodes

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
                        # x direction
                        if con.direction == 1:
                            b[int(x_i), t_n] = self.ramp_loading( deck, t_n , con )
                        # y direction
                        if con.direction == 2:
                            b[int(x_i)+ deck.num_nodes , t_n] = self.ramp_loading( deck, t_n , con )
                        # z direction
                        if con.direction == 3:
                            b[int(x_i) + 2 * deck.num_nodes, t_n] = self.ramp_loading( deck, t_n , con )
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
    def compute_m(self, y,dim,length):
        if dim == 1:
            M = np.zeros((length, length))
            for x_i in range(0, length):
                index_x_family = self.neighbors.get_index_x_family(x_i)
                for x_p in index_x_family:
                    M[x_i, x_p] = (y[x_p] - y[x_i]) / np.absolute(y[x_p] - y[x_i])
            return M
        if dim == 2:
            M_x = np.zeros((length, length))
            M_y = np.zeros((length, length))
            for x_i in range(0, length):
                index_x_family = self.neighbors.get_index_x_family(x_i)
                for x_p in index_x_family:
                    distance = np.sqrt(np.power(y[x_p] - y[x_i],2)+np.power(y[length +x_p] - y[length + x_i],2))
                    M_x[x_i, x_p] = (y[x_p] - y[x_i]) / distance
                    M_y[x_i, x_p] = (y[length + x_p] - y[ length + x_i]) / distance
            return M_x , M_y
        if dim == 3:
            M_x = np.zeros((length, length))
            M_y = np.zeros((length, length))
            M_z = np.zeros((length, length))
            for x_i in range(0, length):
                index_x_family = self.neighbors.get_index_x_family(x_i)
                for x_p in index_x_family:
                    distance = np.sqrt(np.power(y[x_p] - y[x_i],2) + np.power(y[length +x_p] - y[length + x_i],2) + np.power(y[2*length +x_p] - y[2*length + x_i],2))
                    M_x[x_i, x_p] = (y[x_p] - y[x_i]) / distance
                    M_y[x_i, x_p] = (y[length + x_p] - y[ length + x_i]) / distance
                    M_z[x_i, x_p] = (y[2*length + x_p] - y[ 2*length + x_i]) / distance
            return M_x , M_y , M_z
        
    # Computes the weights for each PD node
    def weighted_function(self, deck, x, x_i):
        index_x_family = self.neighbors.get_index_x_family(x_i)
        #print "inside", x_i, index_x_family
        result = 0
        for x_p in index_x_family:
            if deck.dim == 1:
                result += deck.influence_function * (x[x_p] - x[x_i])**2 * deck.geometry.volumes[x_p]
            if deck.dim == 2:
                actual = np.power(x[x_p][0] - x[x_i][0],2)+np.power(x[x_p][1] - x[x_i][1],2)
                result += deck.influence_function * actual * deck.geometry.volumes[x_p]
            if deck.dim == 3:
                actual = np.power(x[x_p][0] - x[x_i][0],2)+np.power(x[x_p][1] - x[x_i][1],2)+np.power(x[x_p][2] - x[x_i][2],2)
                result += deck.influence_function * actual * deck.geometry.volumes[x_p]
        return result

    # Computes the residual vector used in the quasi_static_solver function
    def compute_residual(self, y, deck, t_n):
        residual = np.zeros( ( self.len_x ) )
        for con in deck.conditions:
            if con.type == "Displacement": 
                for id_node in con.id:
                    # x direction
                    if con.direction == 1:
                        if deck.dim == 1:
                            y[int(id_node)] = deck.geometry.nodes[int(id_node)] + con.value
                        else:
                            y[int(id_node)] = deck.geometry.nodes[int(id_node)][0] + con.value
                    # y direction
                    if con.direction == 2:
                        y[int(id_node)+ deck.num_nodes] = deck.geometry.nodes[int(id_node)][1] + con.value    
                    # z direction
                    if con.direction == 3:
                        y[int(id_node)+ 2 * deck.num_nodes] = deck.geometry.nodes[int(id_node)][2] + con.value
                            
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

    # This function solves the problem at each time step, using the previous
    # time step solution as an initial guess.
    # This function calls the compute_residual function
    def quasi_static_solver(self, y, deck):
        
        for t_n in range(1, deck.time_steps):
            solver = scipy.optimize.root(self.compute_residual, y, args=(deck, t_n), method=deck.solver_type,jac=None,tol=deck.solver_tolerance,callback=None,options={'maxiter':1000,'xtol':1.0e-12,'xatol':1.0e-12,'ftol':1.0e-12})
            if deck.dim == 1:
                self.y[:, t_n] = solver.x[:,0]
            else:
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
        if deck.dim == 1:
            y = z + 0.25 * random.uniform(-1, 1) * deck.delta_x
        if deck.dim >= 2:
            y[:deck.num_nodes] = z[:deck.num_nodes] + 0.25 * random.uniform(-1, 1) * deck.delta_x
            y[deck.num_nodes:2*deck.num_nodes] = z[deck.num_nodes:2*deck.num_nodes] + 0.25 * random.uniform(-1, 1) * deck.delta_y
        if deck.dim >= 3:
            y[2*deck.num_nodes:3*deck.num_nodes] = z[2*deck.num_nodes:3*deck.num_nodes] + 0.25 * random.uniform(-1, 1) * deck.delta_z
        #print y
        return y

    def strain_center_bar(self, deck):
        Mid_Node_1 = int(deck.num_nodes/2)-1
        #print "Mid_Node_1 =" , Mid_Node_1
        Mid_Node_2 = int(deck.num_nodes/2)+1
        #print "Mid_Node_2 =" , Mid_Node_2
        for t_n in range(1, deck.time_steps):
            self.strain[t_n] = (np.absolute(self.y[Mid_Node_2,t_n] - self.y[Mid_Node_1,t_n]) - np.absolute(deck.geometry.nodes[Mid_Node_2][0] - deck.geometry.nodes[Mid_Node_1][0])) / np.absolute(deck.geometry.nodes[Mid_Node_2][0] - deck.geometry.nodes[Mid_Node_1][0])
    
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
#            for x_i in range(0, deck.num_nodes):               
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
#            for x_i in range(0, deck.num_nodes):
#                index_x_family = self.get_index_x_family(x_i)
#                modulus = deck.get_elastic_material_properties()
#                for x_p in index_x_family:
#                    #Silling-Askari2005, Eq17
#                    stretch = (abs(self.u[x_p, t_n] - self.u[x_i, t_n]) - abs(
#                        deck.geometry.nodes[x_p] - deck.geometry.nodes[x_i])) / abs(deck.geometry.nodes[x_p] - deck.geometry.nodes[x_i])
#                    #Silling-Askari2005, Eq22
#                    energy[x_i, t_n] = (
#                        math.pi * modulus * math.pow(stretch, 2) * math.pow(self.Horizon, 4)) / 4.
#        self.strain_energy = energy
