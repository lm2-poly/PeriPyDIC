# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 12:50:28 2015

@author: ilyass.tabiai@gmail.com
@author: rolland.delorme@gmail.com
@author: diehl@ins.uni-bonn.de
"""

import logging
from scipy.optimize import fsolve
import timeit
#from deck_elas import PD_deck
import IO.PD_deck
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
import random
import csv
import pdb
import math
import datetime

logger = logging.getLogger(__name__)


class PD_problem():

    def __init__(self, PD_deck):
        # Import initial data
        self.x = np.zeros(PD_deck.Num_Nodes)
        self.y = np.zeros((int(PD_deck.Num_Nodes), int(PD_deck.Num_TimeStep)))
        self.u = np.zeros((int(PD_deck.Num_Nodes), int(PD_deck.Num_TimeStep)))
        self.energy = np.zeros(
            (int(PD_deck.Num_Nodes), int(PD_deck.Num_TimeStep)))
        self.y[:, 0] = self.x
        self.strain = np.zeros((int(PD_deck.Num_TimeStep)))
        self.forces = np.zeros(
            (int(PD_deck.Num_Nodes), int(PD_deck.Num_TimeStep)))
        self.ext = np.zeros(
            (int(
                PD_deck.Num_Nodes), int(
                PD_deck.Num_Nodes), int(
                PD_deck.Num_TimeStep)))
        self.experimental_nodes = 3
        self.exp_displacement = np.zeros(
            (int(PD_deck.Num_TimeStep) - 1, self.experimental_nodes))
        self.exp_times = np.zeros((int(PD_deck.Num_TimeStep) - 1))
        self.exp_init_positions = np.zeros(self.experimental_nodes)

        self.get_pd_nodes(PD_deck)
        if PD_deck.Loading_Flag == "RAMP":
            self.compute_b(PD_deck)
        elif PD_deck.Loading_Flag == "LINEAR_DISPLACEMENT":
            self.compute_u_load(PD_deck)
        self.compute_horizon(PD_deck)
        self.generate_neighborhood_matrix(PD_deck, self.x)

        # For viscoelasticity
        #self.Modulus, self.Relaxation_Time = PD_deck.get_viscoelastic_material_properties()
        #self.ext_visco = np.zeros( ( int(PD_deck.Num_Nodes), int(PD_deck.Num_Nodesodesodes), len(self.Relaxation_Time), int(PD_deck.Num_TimeStep)) )

    # Creates a loading vector b which describes the force applied on each node
    # at any time step
    def compute_b(self, PD_deck):
        # Build  matrix b[row = node, column = time]
        b = np.zeros((int(PD_deck.Num_Nodes), int(PD_deck.Num_TimeStep)))
        PD_deck.get_parameters_loading_ramp()
        if PD_deck.Loading_Flag == "RAMP":
            for x_i in range(len(self.x) -
                             PD_deck.Horizon_Factor, len(self.x)):
                for t_n in range(1, int(PD_deck.Num_TimeStep)):
                    b[x_i, t_n] = self.ramp_loading(PD_deck, t_n)
        else:
            logger.error(
                "There is a problem with the Boundary Conditions in your XML deck.")
        # print b
        self.b = b

    # Adds a displacement step to the points on the edges of the bar
    # Only additional points' (Horizon_Factor) on the edge of the bar are
    # concerned by this function
    def compute_u_load(self, PD_deck):
        # We only update the nodes we actually pull (Horizon_Factor)
        displ_load = np.zeros(
            (int(
                PD_deck.Horizon_Factor), int(
                PD_deck.Num_TimeStep)))
        PD_deck.get_parameters_linear_displacement()
        for x_i in range(0, PD_deck.Horizon_Factor):
            for t_n in range(0, int(PD_deck.Num_TimeStep)):
                
                displ_load[
                    x_i, t_n] = self.linear_displacement_loading(
                    PD_deck, t_n)
        self.displ_load = displ_load

    # Creates a vector of linearly distributed nodes along the bar
    def get_pd_nodes(self, PD_deck):
        # Define x
        for i in range(0, PD_deck.Num_Nodes):
            self.x[i] = i * PD_deck.Delta_x

    # Provides ramp force values to compute the load vector b
    def ramp_loading(self, PD_deck, t_n):
        Time_t = PD_deck.Delta_t * t_n
        if Time_t <= PD_deck.Ramp_Time:
            result = (PD_deck.Force_Density * Time_t) / PD_deck.Ramp_Time
            return result
        else:
            result = PD_deck.Force_Density
            return result

    # Provides linear displacement loading applied on edges of bar
    def linear_displacement_loading(self, PD_deck, t_n):
        Time_t = PD_deck.Delta_t * t_n
        displacement = Time_t * PD_deck.Speed
        return displacement

    # Computes the horizon
    def compute_horizon(self, PD_deck):
        # Be sure that points are IN the horizon
        safety_small_fraction = 1.01
        self.Horizon = PD_deck.Horizon_Factor * \
            PD_deck.Delta_x * safety_small_fraction

    # Returns a list of addresses of the neighbors of a point x_i
    def get_index_x_family(self, x_i):
        return (np.where(self.family[x_i] == 1))[0]

    # Generates matrix neighborhood
    def generate_neighborhood_matrix(self, PD_deck, x):
        self.family = np.zeros((int(
                PD_deck.Num_Nodes), int(
                PD_deck.Num_Nodes)))
        for x_i in range(0, len(x)):
            for x_p in range(0, len(x)):
                if x_p == x_i:
                    pass
                elif np.absolute(x_i - x_p) <= PD_deck.Horizon_Factor:
                    self.family[x_i][x_p] = 1
                else:
                    pass

    # Computes the shape tensor (here a scalar) for each node
    def compute_m(self, Num_Nodes, y):
        M = np.zeros((int(Num_Nodes), int(Num_Nodes)))
        for x_i in range(0, len(self.x)):
            index_x_family = self.get_index_x_family(x_i)
            for x_p in index_x_family:
                M[x_i, x_p] = (y[x_p] - y[x_i]) / np.absolute(y[x_p] - y[x_i])
        return M

    # Computes the weights for each PD node
    def weighted_function(self, PD_deck, x, x_i):
        Horizon = self.compute_horizon(PD_deck)
        index_x_family = self.get_index_x_family(x_i)
        Delta_V = PD_deck.Volume
        Influence_Function = PD_deck.Influence_Function
        result = 0
        for x_p in index_x_family:
            result = result + Influence_Function * \
                (x[x_p] - x[x_i]) ** 2 * Delta_V
        return result



    # Comutes the residual vector used in the quasi_static_solver function
    def compute_residual(self, y, PD_deck, t_n):
        residual = np.zeros((int(PD_deck.Num_Nodes)))
        from elastic import elastic_material
        
        variables = elastic_material(PD_deck, self, y)
        self.update_force_data(variables, t_n)
        self.update_ext_state_data(variables, t_n)

        for x_i in range(PD_deck.Horizon_Factor, len(self.x)):
            if PD_deck.Loading_Flag == "RAMP":
                residual[x_i] = variables.Ts[x_i] + self.b[x_i, t_n]
            else:
                residual[x_i] = variables.Ts[x_i]
                # From johntfoster/1DPDpy, the residual should be set to 0 on
                # Boundary conditions
                for x_i in range(0, PD_deck.Horizon_Factor):
                    #Clamped node
                    residual[x_i] = 0
                    #Pulling nodes
                    residual[len(self.x) - 1 - x_i] = 0
        return residual

    # This functtion solves the problem at each time step, using the previous
    # time step solution as an initial guess
    # This function calls the compute_residual function
    def quasi_static_solver(self, y, PD_deck):

        for t_n in range(1, PD_deck.Num_TimeStep):

            if PD_deck.Loading_Flag == "LINEAR_DISPLACEMENT":
                if t_n == 1:
                    for x_i in range(0, PD_deck.Horizon_Factor):
                        self.y[:, t_n] = self.x

            y = self.y[:, t_n]
            
            
            # Update the boundary conditions
            # Clamped Nodes
            for x_i in range(0, PD_deck.Horizon_Factor):
                y[x_i] = self.x[x_i]
            # Nodes being pulled
            for x_i in range(0, PD_deck.Horizon_Factor):
                y[len(self.x) - 1 - x_i] = self.x[len(self.x) -
                                                  1 - x_i] + self.displ_load[x_i, t_n]
            # Compute the forces and energy before solving the equilibrium
            # Else the force is null since equilibrium has been reached
            from elastic import elastic_material
            variables = elastic_material(PD_deck, self, y)
            self.update_energy_data(variables, t_n)
            
            print "t_n", t_n
            print "BEFORE:", y

            solver = scipy.optimize.root(
                self.compute_residual,
                y,
                args=(
                    PD_deck,
                    t_n),
                method='krylov',
                jac=None,
                tol=1.0e-12,
                callback=None,
                options={
                    'maxiter': 1000,
                    'xtol': 1.0e-12,
                    'xatol': 1.0e-12,
                    'ftol': 1.0e-12})

            print "AFTER:", solver.x

            # pdb.set_trace()

            if t_n + 1 > PD_deck.Num_TimeStep - 1:
                pass
            else:
                self.y[:, t_n + 1] = solver.x

            #y = solver.x + 0.1*random.uniform(-1,1)*PD_deck.Delta_x
            y = self.random_initial_guess(solver.x, PD_deck)

            # Notification if the solver failed
            if solver.success == "False":
                logger.warning("Convergence could not be reached.")
            else:
                logger.info(t_n, solver.success)

            self.update_displacements(t_n)

        return solver



    # Records the force vector at each time step
    def update_force_data(self, variables, t_n):
        self.forces[:, t_n] = variables.Ts

    # Records the force vector at each time step
    def update_energy_data(self, variables, t_n):
        self.energy = self.strain_energy_from_force

    # Records the ext_state vector at each time step
    def update_ext_state_data(self, variables, t_n):
        self.ext[:, :, t_n] = variables.e

    # Records the ext_state vector at each time step
    def update_displacements(self, t_n):
        self.u[:, t_n] = self.y[:, t_n] - self.y[:, 0]

    def random_initial_guess(self, z, PD_deck):
        y = np.zeros((int(PD_deck.Num_Nodes)))
        y = z + 0.1 * random.uniform(-1, 1) * PD_deck.Delta_x
        return y

    # Computes the strain energy density from Ts x u
    def strain_energy_from_force(self, PD_deck):
        energy = np.zeros( (int(PD_deck.Num_TimeStep), int(PD_deck.Num_Nodes)) )
        for t_n in range(0, PD_deck.Num_TimeStep ):   
            for x_i in range(0, PD_deck.Num_Nodes):               
                energy[t_n , x_i] = abs(self.forces[x_i, t_n]) * abs(self.u[x_i, t_n]) * PD_deck.Volume
                #abs(self.forces[x_i, t_n]), abs(self.u[x_i, t_n]), energy[t_n , x_i]
        print "ENERGY:",
        print energy
        self.strain_energy_from_force = energy         
    
    #Computes the strian energy using the formula iven in the PMB
    def strain_energy_bond_based(self, PD_deck):
        energy = np.zeros((int(PD_deck.Num_Nodes), int(PD_deck.Num_TimeStep)))
        for t_n in range(0, PD_deck.Num_TimeStep):
            for x_i in range(0, PD_deck.Num_Nodes):
                index_x_family = self.get_index_x_family(x_i)
                modulus = PD_deck.get_elastic_material_properties()
                for x_p in index_x_family:
                    #Silling-Askari2005, Eq17
                    stretch = (abs(self.u[x_p, t_n] - self.u[x_i, t_n]) - abs(
                        self.x[x_p] - self.x[x_i])) / abs(self.x[x_p] - self.x[x_i])
                    #Silling-Askari2005, Eq22
                    energy[x_i, t_n] = (
                        math.pi * modulus * math.pow(stretch, 2) * math.pow(self.Horizon, 4)) / 4.
        self.strain_energy = energy

    # Exports the data to a CSV file
    # Sill needs work...

    def strain_center_bar(self, PD_deck):
        Mid_Node_1 = int(PD_deck.Num_Nodes / 2) - 1
        Mid_Node_2 = int(PD_deck.Num_Nodes / 2)
        for t_n in range(1, PD_deck.Num_TimeStep):
            self.strain[t_n] = (np.absolute(self.y[Mid_Node_2, t_n] - self.y[Mid_Node_1, t_n]) - np.absolute(
                self.x[Mid_Node_2] - self.x[Mid_Node_1])) / np.absolute(self.x[Mid_Node_2] - self.x[Mid_Node_1])

    def write_data_to_csv(self, PD_deck):

        f = open('data_csv', 'wr+')

        f.write("Time,0,Position")
        for node in self.x:
            f.write("," + str(node))
        f.write("\n")

        for t_n in range(1, PD_deck.Num_TimeStep):

            f.write("Time")
            # for node in self.x:
            f.write("," + str(t_n * PD_deck.Delta_t))
            f.write(",Position")
            for position in self.y[:, t_n]:
                f.write("," + str(position))
            f.write("\n")
        f.write("\n")

        f.write("Time,0,Strain,0.0")
        f.write("\n")

        for t_n in range(1, PD_deck.Num_TimeStep):

            f.write("Time")
            # for node in self.x:
            f.write("," + str(t_n * PD_deck.Delta_t))
            f.write(",Strain")
            f.write("," + str(self.strain[t_n]))
            f.write("\n")

            # f.write("Extension_state")
            # for extension in self.ext[:, :, t_n]:
            #    f.write(","+str(extension))
            # f.write("\n")

            # f.write("Force")
            # for force in self.forces[:, t_n]:
            #    f.write(","+str(force))
            # f.write("\n")
            # f.write("\n")

    def read_csv_from_dic(self, file):
        with open(file, 'r') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')

            next(reader, None)
            i = 0
            for row in reader:
                if float(row[0]) == -1:
                    self.exp_init_positions = map(float, row[1:])
                else:

                    self.exp_displacement[i] = np.array(map(float, row[1:]))
                    self.exp_times[i] = float(row[0])
                    i += 1

    def plot_force(self, PD_deck):
        for t_n in range(1, PD_deck.Num_TimeStep):
            force_plot = plt
            force_plot.plot(self.y[:, 1], self.forces[:, t_n], '-+', label=t_n)
        force_plot.legend(title="forces")
        return force_plot

    def plot_positions(self, PD_deck):
        for t_n in range(1, PD_deck.Num_TimeStep):
            position_plot = plt
            position_plot.plot(self.y[:, 1], self.y[:, t_n], '-+')
        position_plot.legend(title="position")
        return position_plot
    
    #For simulaiton, we add Horizon_Factor points to the right of the bar to
    #pull and to the left which are fixed. They should be removed before 
    #presenting the results
    def remove_additional_points(self, PD_deck, energy_list):
        cleaned_energy_list = []
        for t_n in range(0, len(energy_list)):
            energy_t_n = []
            for x_i in range(PD_deck.Horizon_Factor, len(energy_list[t_n])-PD_deck.Horizon_Factor):
                energy_t_n.append(energy_list[t_n][x_i])
            cleaned_energy_list.append(energy_t_n)
        return cleaned_energy_list
            
        

        
    def plot_energy(self,energy,time,initial,outpath):
        #print len(energy) , len(time), len(initial) 
        maxvalues = []
        color = []
        for i in range(0, len(initial)):
            e = []
            for j in range(0, len(energy)):
                e.append(abs(energy[j][i]))
            line = plt.plot(time, e, marker="o")
            color.append(line[0].get_color())
            maxvalues.append(max(e))
        #print color
        #plt.plot([0, max(time)], [max(maxvalues) + 100,
        #                         max(maxvalues) + 100], lw=2, c='black')
        #for i in range(len(initial)):
        #    plt.plot(
        #        (max(time) / 37.5) * initial[i],
        #        max(maxvalues) + 100,
        #        color[i],
        #        marker='o')
        plt.grid()
        plt.xlabel("Time [s]")
        plt.ylabel("Strain Energy")
        plt.savefig(outpath +
                    "energy_" +
                    str(datetime.datetime.now().strftime("%Y-%m-%d_%H:%M")) +
                    ".pdf")
