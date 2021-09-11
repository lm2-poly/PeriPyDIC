# -*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrickdiehl@lsu.edu
import sys
import numpy as np
from scipy import linalg
np.set_printoptions(precision=8, threshold=sys.maxsize)
from ..util import functions

## Class to compute the well-known strain and stress tensors defined in the classical continuum mechanics
class CCM_calcul():

    ## Constructor
    # @param deck The input deck
    # @param data_solver Data from the peridynamic problem/solving class
    def __init__(self, deck, data_solver):

        ## Nodes' initial position
        self.x = deck.geometry.nodes

        ## Nodes' positions stored for each time step
        self.y = data_solver.y

        ## Global internal force density array storing the force density attached to each node for each time step
        self.force_int = data_solver.force_int

        ## Extension array storing the extension for each node between itself and its family
        self.ext = data_solver.ext

        ## Dimension of the data_solver (1D, 2D or 3D)
        self.dim = deck.dim

        ## Amount of nodes in the data_solver
        self.num_nodes = deck.num_nodes

        ## Amount of time step
        self.time_steps = deck.time_steps

        ## Volume related to each node
        self.node_volumes = deck.geometry.volumes

        ## Weighted volume
        self.Weighted_Volume = data_solver.weighted_volume

        ## Volume correction factor
        self.Volume_Correction = data_solver.volume_correction

        ## Material type
        self.material_type = deck.material_type

        if self.material_type == "Elastic":
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

        if self.material_type == "Viscoelastic":
            if deck.dim >= 2:
                ## Bulk modulus of the material
                self.K = deck.relax_bulk_modulus
                ## Shear modulus of the material
                self.Mu = deck.relax_shear_modulus

                if deck.dim == 2:
                    ## Poisson ratio of the material
                    self.Nu = (3. * self.K - 2. * self.Mu) / (2. * (3. * self.K + self.Mu))
                    if deck.type2d == "Plane_Stress":
                        ## Factor applied for 2D plane stress to compute dilatation and force state
                        self.factor2d = (2. * self.Nu - 1.) / (self.Nu - 1.)
                    if deck.type2d == "Plane_Strain":
                        ## Plane strain
                        self.factor2d = 1

        ## Compute the global strain tensor storing the strain tensor for each node at each time step
        self.compute_global_strain_tensor(deck,data_solver)

        ## Compute the displacement for each node at each time step
        #self.compute_u_displacement()

        ## Compute the global stress tensor storing the strain tensor for each node at each time step
        if self.material_type == "Elastic":
            self.compute_global_stress_tensor(deck,data_solver)

    ## Provide the image of (xi - xp) under the reference position vector state X
    # @param data_solver Data from the peridynamic problem/solving class
    # @param i Id of Node "i"
    # @param p Id of Node "p" with Node "i" family
    # @return Image of (xi - xp) under the deformation vector state X
    def X_vector_state(self, data_solver, i, p):
        X = self.x[p,:] - self.x[i,:]
        X = np.reshape(X,(self.dim,1))
        return X

    ## Provide the image of (xi - xp) under the deformation vector state Y
    # @param data_solver Data from the peridynamic problem/solving class
    # @param i Id of Node "i"
    # @param p Id of Node "p" with Node "i" family
    # @param t_n Id of the time step
    # @return Image of (xi - xp) under the deformation vector state Y
    def Y_vector_state(self, data_solver, i, p, t_n):
        Y = self.y[p,:,t_n] - self.y[i,:,t_n]
        Y = np.reshape(Y,(self.dim,1))
        return Y

    ## Provide the shape tensor K related to Node "i"
    # @param data_solver Data from the peridynamic problem/solving class
    # @param deck Input deck
    # @param i Id of Node "i"
    # @return Shape tensor K
    def K_shape_tensor(self, deck, data_solver, i):
        K = np.zeros((self.dim, self.dim),dtype=np.float64)
        index_x_family = data_solver.neighbors.get_index_x_family(i)
        n = 0
        for p in index_x_family:
            X = self.X_vector_state(data_solver, i, p)
            K += functions.w(data_solver, X, deck.influence_function) * np.dot(X,X.T) * self.Volume_Correction[i,n] * self.node_volumes[p]
            n += 1
        return K

    ## Provide the deformation gradient tensor related to Node "i"
    # @param deck Input deck
    # @param data_solver Data from the peridynamic problem/solving class
    # @param i Id of Node "i"
    # @param t_n Id of the time step
    # @return Deformation gradient tensor related to Node "i"
    def deformation_gradient(self, deck, data_solver, i, t_n):
        tmp = np.zeros((self.dim, self.dim),dtype=np.float64)
        index_x_family = data_solver.neighbors.get_index_x_family(i)
        n = 0
        for p in index_x_family:
            Y = self.Y_vector_state(data_solver, i, p, t_n)
            X = self.X_vector_state(data_solver, i, p)
            tmp += functions.w(data_solver, X, deck.influence_function) * np.dot(Y,X.T) * self.Volume_Correction[i,n] * self.node_volumes[p]
            n += 1
        deformation = np.dot(tmp, linalg.inv(self.K_shape_tensor(deck,data_solver, i)))
        return deformation

    ## Provide the strain tensor related to Node "i"
    # @param data_solver Data from the peridynamic problem/solving class
    # @param i Id of Node "i"
    # @param t_n Id of the time step
    # @return strain tensor related do Node "i"
    def strain_tensor(self, deck, data_solver, i, t_n):
        F = self.deformation_gradient(deck,data_solver, i, t_n)
        strain = (F + F.T)/2 - np.identity(self.dim, dtype=np.float64)
        return strain

    ## Compute the global strain tensor storing the strain tensor for each node at each time step
    # @param data_solver Data from the peridynamic problem/solving class
    def compute_global_strain_tensor(self, deck, data_solver):
        ## Golbal strain tensor storing the strain tensor for each node at each time step
        self.global_strain = np.zeros((self.num_nodes*self.dim, self.dim, self.time_steps),dtype=np.float64)
        for t_n in range(1, self.time_steps):
            for i in range(0, self.num_nodes):
                tmp = self.strain_tensor(deck,data_solver, i, t_n)
                for j in range(0, self.dim):
                    for r in range(0, self.dim):
                        self.global_strain[i*self.dim+r,j,t_n] = tmp[r,j]

    ## Provide the image of x under the Dirac Delta Function
    # @param x Vector x
    # @return 1 if x is a null-vector, otherwise 0
    def DiracDelta(self, x, i, q, m):
        if linalg.norm(x) == 0.:
            delta = 1. / (self.node_volumes[q] * self.Volume_Correction[i,m])
        else:
            delta = 0.
        return delta

    ## Provide the modulus state K related to Node "i"
    # @param data_solver Data from the peridynamic problem/solving class
    # @param i Id of Node "i"
    # @param p Id of Node "p" with Node "i" family
    # @param q Id of Node "q" with Node "i" family
    # @param t_n Id of the time step
    # @return Shape tensor K
    def K_modulus_tensor(self, deck, data_solver, i , p, q, m):
        Xp = self.X_vector_state(data_solver, i, p)
        M = Xp / linalg.norm(Xp)
        Xq = self.X_vector_state(data_solver, i, q)
        xp = deck.geometry.nodes[p,:] - deck.geometry.nodes[i,:]
        xq = deck.geometry.nodes[q,:] - deck.geometry.nodes[i,:]
        if self.material_type == "Elastic":
            if self.dim == 1:
                # PD material parameter
                alpha = self.Young_Modulus / self.Weighted_Volume[i]
                K = alpha * functions.w(data_solver, xp, deck.influence_function) * np.dot(M,M.T) * self.DiracDelta(Xq - Xp, i, q, m)

            if self.dim == 2:
                # PD material parameter
                # Plane stress
                alpha_s = (9. / self.Weighted_Volume[i]) * (self.K + ((self.Nu + 1.)/(2. * self.Nu - 1.))**2 * self.Mu / 9.)
                # Plane strain
                #alpha_s = (9. / self.Weighted_Volume[i]) * (self.K + self.Mu / 9.)
                alpha_d = (8. / self.Weighted_Volume[i]) * self.Mu
                alpha_sb = (2. * self.factor2d * alpha_s - (3. - 2. * self.factor2d) * alpha_d) /3.
                K = ((alpha_sb - alpha_d) / self.Weighted_Volume[i]) * functions.w(data_solver, xp, deck.influence_function) * functions.w(data_solver, xq, deck.influence_function) * np.dot(Xp,Xq.T) + alpha_d * functions.w(data_solver, xp, deck.influence_function) * np.dot(M,M.T) * self.DiracDelta(Xq - Xp, i, q, m)

            if self.dim == 3:
                # PD material parameter
                alpha_s = (9. / self.Weighted_Volume[i]) * self.K
                alpha_d = (15. / self.Weighted_Volume[i]) * self.Mu
                K = ((alpha_s - alpha_d) / self.Weighted_Volume[i]) * functions.w(data_solver, xp, deck.influence_function) * functions.w(data_solver, xq, deck.influence_function) * np.dot(Xp,Xq.T) + alpha_d * functions.w(data_solver, xp, deck.influence_function) * np.dot(M,M.T) * self.DiracDelta(Xq - Xp, i, q, m)

        return K

    ## Provide the stress tensor related to Node "i"
    # @param data_solver Data from the peridynamic problem/solving class
    # @param i Id of Node "i"
    # @param t_n Id of the time step
    # @return stress tensor related do Node "i"
    def stress_tensor(self, deck, data_solver, i, t_n):
        #force = np.zeros((self.dim, 1),dtype=np.float64)
        stress = np.zeros((self.dim, self.dim),dtype=np.float64)
        index_x_family = data_solver.neighbors.get_index_x_family(i)
        n = 0
        for p in index_x_family:
            Xp = self.X_vector_state(data_solver, i, p)
            m = 0
            for q in index_x_family:
                Xq = self.X_vector_state(data_solver, i, q)
                stress += np.dot(np.dot(self.K_modulus_tensor(deck,data_solver, i , p, q, m), np.dot(self.strain_tensor(deck,data_solver, i, t_n), Xq)), Xp.T) * self.Volume_Correction[i,m] * self.node_volumes[q] * self.Volume_Correction[i,n] * self.node_volumes[p]
                m += 1
            n += 1
        return stress

    ## Compute the global stress tensor storing the strain tensor for each node at each time step
    # @param data_solver Data from the peridynamic problem/solving class
    def compute_global_stress_tensor(self, deck, data_solver):
        ## Golbal strain tensor storing the strain tensor for each node at each time step
        self.global_stress = np.zeros((self.num_nodes*self.dim, self.dim, self.time_steps),dtype=np.float64)
        for t_n in range(1, self.time_steps):
            for i in range(0, self.num_nodes):
                tmp = self.stress_tensor(deck,data_solver, i, t_n)
                for j in range(0, self.dim):
                    for r in range(0, self.dim):
                        self.global_stress[i*self.dim+r,j,t_n] = tmp[r,j]

    ## Compute the displacement for each node at each time step
    def compute_u_displacement(self):
        ## Displacement vector between two consecutives time steps for each node
        self.u = np.zeros((self.num_nodes, self.dim, self.time_steps),dtype=np.float64)
        for t_n in range(1, self.time_steps):
            for i in range(0, self.num_nodes):
                self.u[i,:,t_n] = self.y[i,:,t_n] - self.y[i,:,t_n-1]

    ## Provide the image of (xi - xp) under the displacement vector state U
    # @param i Id of Node "i"
    # @param p Id of Node "p" with Node "i" family
    # @param t_n Id of the time step
    # @return Image of (xi - xp) under the deformation vector state U
    def U_vector_state(self, i, p, t_n):
        U = self.u[p,:,t_n] - self.u[i,:,t_n]
        U = np.reshape(U,(self.dim,1))
        return U
