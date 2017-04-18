# -*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrick.diehl@polymtl.ca
import numpy as np
from scipy import linalg
np.set_printoptions(precision=8, threshold='nan', suppress=True)
from multiprocessing import Process, Lock
import sharedmem
import util.linalgebra
import sys

## Class to compute the global internal volumic force at each node of an elastic material using its material properties
class Viscoelastic_material():

    ## Constructor
    # @param deck The input deck
    # @param data_solver Data from the peridynamic problem/solving class
    # @param y Actual nodes' position
    # @param t_n Id of the time step
    def __init__(self, deck, data_solver, y, t_n):
        
        ## Influence function
        self.w = deck.influence_function

        ## Weighted volume
        self.Weighted_Volume = data_solver.weighted_volume
        
        if deck.dim == 1:
            ## Relaxation modulus of the material
            self.Relax_Modulus = deck.relax_modulus
            ## Relaxation time of the material
            self.Relax_Time = deck.relax_time
        
        if deck.dim == 2:
            print "Error: 2D problem not implemented in viscoelasticity"
            sys.exit(1)

        if deck.dim == 3:
            print "Error: 3D problem not implemented in viscoelasticity"
            sys.exit(1)

        ## Compute the dilatation for each node
        self.compute_dilatation(deck, data_solver, y)
        
        ## Compute the global internal force density at each node
        self.compute_f_int(deck, data_solver, y)
        
    ## Compute the dilatation for each node
    # @param deck The input deck
    # @param data_solver Data from the peridynamic problem/solving class
    # @param y The actual nodes' position
    # @param start Starting Id of the loop
    # @param end Ending Id of the loop
    def compute_dilatation_slice(self, deck, data_solver, y, start, end):
        for i in range(start, end):
            index_x_family = data_solver.neighbors.get_index_x_family(i)
            for p in index_x_family:
                    Y = (y[p,:]) - y[i,:]
                    X = deck.geometry.nodes[p,:] - deck.geometry.nodes[i,:]
                    self.e[i,p] = util.linalgebra.norm(Y) - util.linalgebra.norm(X)
                    
                    if deck.dim == 1:
                        self.dilatation[i] += (1. / self.Weighted_Volume[i]) * self.w * util.linalgebra.norm(X) * self.e[i,p] * deck.geometry.volumes[p]
        
                    if deck.dim == 2:
                        self.dilatation[i] += (2. / self.Weighted_Volume[i]) * self.factor2d * self.w * util.linalgebra.norm(X) * self.e[i,p] * deck.geometry.volumes[p]
        
                    if deck.dim == 3:
                        self.dilatation[i] += (3. / self.Weighted_Volume[i]) * self.w * util.linalgebra.norm(X) * self.e[i,p] * deck.geometry.volumes[p]
                        
    ## Compute the dilatation and also the scalar extension state for each node
    # @param deck The input deck
    # @param data_solver Data from the peridynamic problem/solving class
    # @param y The actual nodes' position                        
    def compute_dilatation(self, deck, data_solver, y):
        ## Dilatation at each node        
        self.dilatation = sharedmem.empty((deck.num_nodes),dtype=np.float64)
        ## Extension between Node "i" and Node "p" within its family
        self.e = sharedmem.empty((deck.num_nodes, deck.num_nodes),dtype=np.float64)  
        
        threads = deck.num_threads
        part = int(deck.num_nodes/threads)
        
        processes = []
        
        for i in range(0,threads):
            start = i * part
            if i < threads - 1:
                end = (i+1) * part
            else:
                end = deck.num_nodes

            processes.append(Process(target=self.compute_dilatation_slice, args=(deck, data_solver, y, start, end)))
            processes[i].start()
           
        for p in processes:
            p.join() 

    ## Compute the viscoelastic part of the scalar extension state
    # @param deck The input deck
    # @param data_solver Data from the peridynamic problem/solving class
    # @param y The actual nodes' position 
    # @param t_n Id of the time step
    # @param start Starting Id of the loop
    # @param end Ending Id of the loop 
    def compute_ext_state_visco_slice(self, deck, data_solver, y, t_n, start, end):
        for i in range(start, end):
            index_x_family = data_solver.neighbors.get_index_x_family(i)
            for p in index_x_family:
                for k in range(1, len(self.Relax_Time)):
                    temp_exp = np.exp((- deck.delta_t)/(self.Relax_Time[k]))
                    delta_e = self.e[x_i, x_p] - problem.ext[x_i, x_p, t_n-1]
                    beta = 1-(self.Relax_Time[k]*(1-temp_exp)) / deck.delta_t
                    e_visco[x_i, x_p, k] = problem.ext[x_i, x_p, t_n-1] *(1-temp_exp) + problem.ext_visco[x_i, x_p, k, t_n-1]*temp_exp + beta*delta_e

    ## Compute the viscoelastic part of the scalar extension state
    # @param deck The input deck
    # @param data_solver Data from the peridynamic problem/solving class
    # @param y The actual nodes' position 
    # @param t_n Id of the time step
    def compute_ext_state_visco(self, deck, data_solver, y, t_n):
        self.e_visco = sharedmem.empty((deck.num_nodes, deck.num_nodes, len(self.Relax_Time)),dtype=np.float64)
        
        threads = deck.num_threads
        part = int(deck.num_nodes/threads)
        
        processes = []
        
        for i in range(0,threads):
            start = i * part
            if i < threads - 1:
                end = (i+1) * part
            else:
                end = deck.num_nodes

            processes.append(Process(target=self.compute_ext_state_visco_slice, args=(deck, data_solver, y, t_n, start, end)))
            processes[i].start()
           
        for p in processes:
            p.join() 



    ## Function to compute the viscous part of the scalar force state
    # @param deck The input deck
    # @param problem The related peridynamic problem
    # @param y The actual postions
    def compute_t_visco(self, deck, problem, y):
        t_visco = np.zeros( (self.len_x, self.len_x ) )
        for x_i in range(0, self.len_x):
            index_x_family = problem.neighbors.get_index_x_family(x_i)
            for x_p in index_x_family:
                for k in range(1, len(self.Relax_Time)):
                    t_visco[x_i, x_p] = t_visco[x_i, x_p] + (self.w / problem.weighted_function(deck, deck.geometry.nodes, x_i))*self.Relax_Modulus[k]*(self.e[x_i, x_p] - self.e_visco[x_i, x_p, k])
        ## Viscous part of the scalar state
        self.t_visco = t_visco

    ## Function to compute the vector force state
    # @param deck The input deck
    # @param problem The related peridynamic problem
    # @param y The actual postions
    def compute_T(self, deck, problem, y):
        tscal = np.zeros( (self.len_x, self.len_x ) )
        for x_i in range(0, self.len_x):
            index_x_family = problem.neighbors.get_index_x_family(x_i)
            for x_p in index_x_family:
                tscal[x_i, x_p] = (self.w / problem.weighted_function(deck, deck.geometry.nodes, x_i))*self.Relax_Modulus[0]*self.e[x_i, x_p] + self.t_visco[x_i, x_p]
        ## Scalar force state
        self.tscal = tscal

        T = np.zeros( (self.len_x, self.len_x ) )
        for x_i in range(0, self.len_x):
            index_x_family = problem.neighbors.get_index_x_family(x_i)
            for x_p in index_x_family:
                T[x_i, x_p] = tscal[x_i, x_p] * self.M[x_i, x_p]
        ## Vector force state
        self.T = T

    ## Function to compute, as a vector state, the global internal volumic force within the equation of motion
    # @param deck The input deck
    # @param problem The related peridynamic problem
    def compute_Ts(self, deck, problem):
        Ts = np.zeros( (self.len_x ) )
        for x_i in range(0, self.len_x):
            index_x_family = problem.neighbors.get_index_x_family(x_i)
            for x_p in index_x_family:
                Ts[x_i] = Ts[x_i] + self.T[x_i, x_p] - self.T[x_p, x_i]
            Ts[x_i] = Ts[x_i] * deck.geometry.volumes[x_i]
        ## Sum of the (T[xi] - T[xp])*Vol[xi] equivalent to the global internal volumic force
        self.Ts = Ts
