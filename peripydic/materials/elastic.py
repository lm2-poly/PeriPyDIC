# -*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrick.diehl@polymtl.ca
import numpy as np
from scipy import linalg
from multiprocessing import Process, Lock
import sharedmem
from ..util import linalgebra
from ..util import functions


## Class to compute the global internal volumic force at each node of an elastic material using its material properties
class Elastic_material():

    ## Constructor
    # @param deck The input deck
    # @param data_solver Data from the peridynamic problem/solving class
    # @param y Actual nodes' position
    def __init__(self, deck, data_solver, y):

        ## Weighted volume
        self.Weighted_Volume = data_solver.weighted_volume

        ## Volume correction factor
        self.Volume_Correction = data_solver.volume_correction

        if deck.dim == 1:
            ## Young modulus of the material
            self.Young_Modulus = deck.young_modulus

        if deck.dim >= 2:
            ## Bulk modulus of the material
            self.K = deck.bulk_modulus
            ## Shear modulus of the material
            self.Mu = deck.shear_modulus

            if deck.dim == 2:
                ## Poisson ratio of the material
                self.Nu = (3. * self.K - 2. * self.Mu) / (2. * (3. * self.K + self.Mu))
                if deck.type2d == "Plane_Stress":
                    ## Factor applied for 2D plane stress to compute dilatation and force state
                    self.factor2d = (2. * self.Nu - 1.) / (self.Nu - 1.)
                if deck.type2d == "Plane_Strain":
                    ## Plane strain
                    self.factor2d = 1

        ## Compute the dilatation for each node
        self.compute_dilatation(deck, data_solver, y)

        ## Compute the global internal force density at each node
        self.compute_f_int(deck, data_solver, y)

        ## Compute the strain energy density at each node
        self.compute_strain_energy(deck, data_solver)

    ## Compute the dilatation for each node
    # @param deck The input deck
    # @param data_solver Data from the peridynamic problem/solving class
    # @param y The actual nodes' position
    # @param start Starting Id of the loop
    # @param end Ending Id of the loop
    def compute_dilatation_slice(self, deck, data_solver, y, start, end):
        for i in range(start, end):
            index_x_family = data_solver.neighbors.get_index_x_family(i)
            n = 0
            for p in index_x_family:
                Y = (y[p,:]) - y[i,:]
                X = deck.geometry.nodes[p,:] - deck.geometry.nodes[i,:]
                self.e[i,n] = linalgebra.norm(Y) - linalgebra.norm(X)

                if deck.dim == 1:
                    self.dilatation[i] += (1. / self.Weighted_Volume[i]) * functions.w(data_solver, X, deck.influence_function) * linalgebra.norm(X) * self.e[i,n] * self.Volume_Correction[i,n] * deck.geometry.volumes[p]

                if deck.dim == 2:
                    if self.Weighted_Volume[i] == 0:
                        print(self.Weighted_Volume[i],i,p)
                    self.dilatation[i] += (2. / self.Weighted_Volume[i]) * self.factor2d * functions.w(data_solver, X, deck.influence_function) * linalgebra.norm(X) * self.e[i,n] * self.Volume_Correction[i,n] * deck.geometry.volumes[p]

                if deck.dim == 3:
                    self.dilatation[i] += (3. / self.Weighted_Volume[i]) * functions.w(data_solver, X, deck.influence_function) * linalgebra.norm(X) * self.e[i,n] * self.Volume_Correction[i,n] * deck.geometry.volumes[p]
                n += 1

    ## Compute the dilatation and and also the scalar extension state for each node
    # @param deck The input deck
    # @param data_solver Data from the peridynamic problem/solving class
    # @param y The actual nodes' position
    def compute_dilatation(self, deck, data_solver, y):
        ## Dilatation at each node
        self.dilatation = sharedmem.empty((deck.num_nodes),dtype=np.float64)
        ## Extension between Node "i" and Node "p" within its family
        self.e = sharedmem.empty((deck.num_nodes, data_solver.neighbors.max_neighbors),dtype=np.float64)

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

    ## Compute the global internal force density at each node
    # @param deck The input deck
    # @param data_solver Data from the peridynamic problem/solving class
    # @param y The actual nodes' position
    # @param start Starting Id of the loop
    # @param end Ending Id of the loop
    def compute_f_int_slice(self, deck, data_solver, y, start, end, data):
        for i in range(start, end):
            index_x_family = data_solver.neighbors.get_index_x_family(i)
            n = 0
            for p in index_x_family:
                Y = y[p,:] - y[i,:]
                X = deck.geometry.nodes[p,:] - deck.geometry.nodes[i,:]

                # Compute the direction vector between Node_p and Node_i
                M = Y / linalgebra.norm(Y)

                if deck.dim == 1:
                    # PD material parameter
                    alpha = self.Young_Modulus / self.Weighted_Volume[i]
                    ## Scalar force state
                    self.t = alpha * functions.w(data_solver, X, deck.influence_function) * self.e[i,n]

                if deck.dim == 2:
                    # PD material parameter
                    if deck.type2d == "Plane_Stress":
                        alpha_s = (9. / self.Weighted_Volume[i]) * (self.K + ((self.Nu + 1.)/(2. * self.Nu - 1.))**2 * self.Mu / 9.)
                    if deck.type2d == "Plane_Strain":
                        alpha_s = (9. / self.Weighted_Volume[i]) * (self.K + self.Mu / 9.)

                    alpha_d = (8. / self.Weighted_Volume[i]) * self.Mu
                    # Scalar extension states
                    e_s = self.dilatation[i] * linalgebra.norm(X) / 3.
                    e_d = self.e[i,n] - e_s
                    # Scalar force states
                    t_s = (2. * self.factor2d * alpha_s - (3. - 2. * self.factor2d) * alpha_d) * functions.w(data_solver, X, deck.influence_function) * e_s / 3.
                    t_d = alpha_d * functions.w(data_solver, X, deck.influence_function) * e_d
                    self.t = t_s + t_d

                if deck.dim == 3:
                    # PD material parameter
                    alpha_s = (9. / self.Weighted_Volume[i]) * self.K
                    alpha_d = (15. / self.Weighted_Volume[i]) * self.Mu
                    # Scalar extension states
                    e_s = self.dilatation[i] * linalgebra.norm(X) / 3.
                    e_d = self.e[i,n] - e_s
                    # Scalar force states
                    t_s = alpha_s * functions.w(data_solver, X, deck.influence_function) * e_s
                    t_d = alpha_d * functions.w(data_solver, X, deck.influence_function) * e_d
                    self.t = t_s + t_d

                #lock.acquire()
                data[i,:] += self.t * M * self.Volume_Correction[i,n] * deck.geometry.volumes[p]
                data[p,:] += -self.t * M * self.Volume_Correction[i,n] * deck.geometry.volumes[i]
                #lock.release()
                n += 1

    ## Compute the global internal force density at each node
    # @param deck The input deck
    # @param data_solver Data from the peridynamic problem/solving class
    # @param y The actual nodes' position
    def compute_f_int(self, deck, data_solver, y):
        ## Internal force density at each node
        self.f_int = sharedmem.empty((deck.num_nodes, deck.dim),dtype=np.float64)
        #self.f_int.fill(0.0)
        #lock = Lock()
        threads = deck.num_threads
        part = int(deck.num_nodes/threads)

        processes = []
        data = []
        for i in range(0,threads):
            start = i * part
            if i < threads - 1:
                end = (i+1) * part
            else:
                end = deck.num_nodes
            data.append(sharedmem.empty((deck.num_nodes, deck.dim),dtype=np.float64))
            #data[i].fill(0)
            processes.append(Process(target=self.compute_f_int_slice, args=(deck, data_solver, y, start, end, data[i])))
            processes[i].start()

        for p in processes:
            p.join()

        for i in range(0,threads):
            self.f_int += data[i]

    ## Computes the strain energy density for each PD node
    # @param deck The input deck
    # @param data_solver Data from the peridynamic problem/solving class
    # @param start Starting Id of the loop
    # @param end Ending Id of the loop
    def compute_strain_energy_slice(self, deck, data_solver, start, end):
        for i in range(start, end):
            index_x_family = data_solver.neighbors.get_index_x_family(i)
            n = 0
            for p in index_x_family:
                X = deck.geometry.nodes[p,:] - deck.geometry.nodes[i,:]
                if deck.dim == 1:
                    # PD material parameter
                    alpha = self.Young_Modulus / self.Weighted_Volume[i]
                    # Strain energy density
                    self.strain_energy[i] += 0.5 * alpha * functions.w(data_solver, X, deck.influence_function) * self.e[i,n]**2 * self.Volume_Correction[i,n] * deck.geometry.volumes[p]

                if deck.dim >= 2:
                    X = deck.geometry.nodes[p,:] - deck.geometry.nodes[i,:]
                    if deck.dim == 2:
                        # PD material parameter
                        if deck.type2d == "Plane_Stress":
                            alpha_s = (9. / self.Weighted_Volume[i]) * (self.K + ((self.Nu + 1.)/(2. * self.Nu - 1.))**2 * self.Mu / 9.)
                        if deck.type2d == "Plane_Strain":
                            alpha_s = (9. / self.Weighted_Volume[i]) * (self.K + self.Mu / 9.)
                        alpha_d = (8. / self.Weighted_Volume[i]) * self.Mu

                    if deck.dim == 3:
                        # PD material parameter
                        alpha_s = (9. / self.Weighted_Volume[i]) * self.K
                        alpha_d = (15. / self.Weighted_Volume[i]) * self.Mu

                    # Scalar extension states
                    e_s = self.dilatation[i] * linalgebra.norm(X) / 3.
                    e_d = self.e[i,n] - e_s
                    # Strain energy density
                    self.strain_energy[i] += 0.5 * functions.w(data_solver, X, deck.influence_function) * (alpha_s * e_s**2 + alpha_d * e_d**2) * self.Volume_Correction[i,n] * deck.geometry.volumes[p]
                n += 1

    ## Compute the strain energy density at each node
    # @param deck The input deck
    # @param data_solver Data from the peridynamic problem/solving class
    def compute_strain_energy(self, deck, data_solver):
        ## Strain energy density at each node
        self.strain_energy = sharedmem.empty((deck.num_nodes),dtype=np.float64)

        threads = deck.num_threads
        part = int(deck.num_nodes/threads)

        processes = []

        for i in range(0,threads):
            start = i * part
            if i < threads - 1:
                end = (i+1) * part
            else:
                end = deck.num_nodes
            processes.append(Process(target=self.compute_strain_energy_slice, args=(deck, data_solver, start, end)))
            processes[i].start()

        for p in processes:
            p.join()
