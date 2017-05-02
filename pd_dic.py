#-*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrick.diehl@polymtl.ca
import sys
import getopt
import numpy as np
np.set_printoptions(precision=8, threshold='nan', suppress=True)
import time
from peripydic import *



def main(argv):
    """
    Main
    """
    helptext = sys.argv[0] + " -i input.yaml -t type"
    types = ['pd', 'dic']

    if len(sys.argv) != 5:
        print helptext
        sys.exit(1)

    try:
        opts, args = getopt.getopt(
            argv, "hi:o:t:", ["ifile=","type="])
    except getopt.GetoptError:
        print helptext
        sys.exit(0)

    for opt, arg in opts:
        if opt in ("-i", "--ifile"):
            inputFile = arg
        elif opt in ("-t", "--type"):
            typeIn = arg
    if typeIn not in types:
        print("Error: Only pd or dic types are supported")
        sys.exit(1)

    if typeIn == types[0]:
        deck = IO.deck.PD_deck(inputFile)
        if deck.material_type == "Elastic":
            simulation(deck)
        elif deck.material_type == "Viscoelastic":
            simulation(deck)
        else:
            print "Error in pd_dict.py: Material type unknown, please use Elastic or Viscoelastic"

    if typeIn == types[1]:
        deck = IO.deck.DIC_deck(inputFile)
        dic(deck)

def dic(deck):
    pb_solver_class = DIC_problem(deck)
    ccm_class = IO.ccm.CCM_calcul(deck, pb_solver_class)

    if deck.vtk_writer.vtk_enabled == True:
        deck.vtk_writer.write_data(deck,pb_solver_class,ccm_class)

def simulation(deck):
    t0 = time.time()
    y_0 = deck.geometry.nodes.copy()
    pb_solver_class = problem.pd.PD_problem(deck)
    pb_solver_class.quasi_static_solver(deck, y_0)
    ccm_class = IO.ccm.CCM_calcul(deck, pb_solver_class)

    writeCSV(deck,pb_solver_class)
    if deck.vtk_writer.vtk_enabled == True:
        deck.vtk_writer.write_data(deck,pb_solver_class,ccm_class)

    print "delta_x =" , deck.delta_X
    print "Horizon =" , pb_solver_class.neighbors.horizon

    strain_tensor = ccm_class.global_strain[:,:,deck.time_steps-1]
    print "epsilon_tensor"
    print strain_tensor

    if deck.material_type == "Elastic":
        stress_tensor = ccm_class.global_stress[:,:,deck.time_steps-1]
        print "stress_tensor"
        print stress_tensor

    strain_longi = pb_solver_class.strain_calculation(deck, 5, 7)
    print "strain_longi", strain_longi
    #print "Nodes positions = "
    #print pb_solver_class.y
    print "Duration:", (time.time() - t0)/60. , "minutes"

def writeCSV(deck,problem):
    for out in deck.outputs:
        if out.outType == "CSV":
            out.write(deck,problem)


# Start the function __main__ at __init__ call
if __name__ == "__main__":
    main(sys.argv[1:])
