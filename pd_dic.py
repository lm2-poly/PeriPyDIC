#-*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrick.diehl@polymtl.ca
import sys
import getopt
import IO.deck
import problem
import numpy as np
np.set_printoptions(threshold='nan')
import time
import pdb

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

def simulation(deck):
    t0 = time.time()
    pb_class = problem.PD_problem(deck)
    y_0 = deck.geometry.nodes.copy()
    pb_class.quasi_static_solver(y_0, deck)
    eps_longi = pb_class.strain_calculation( 3, 5, deck )
    eps_trans = pb_class.strain_calculation( 12, 20, deck )

    writeCSV(deck,pb_class)
    if deck.vtk_writer.vtk_enabled == True:
       deck.vtk_writer.write_data(deck,pb_class)

    print "delta_x =" , deck.delta_X
    print "Horizon =" , pb_class.neighbors.horizon
    print "Strain Longi = " , np.around(eps_longi,decimals=6)
    print "Strain Trans = " , np.around(eps_trans,decimals=6)
    #print "Nodes positions = "
    #print pb_class.y
    print "Duration:", time.time() - t0 , "seconds"

def writeCSV(deck,problem):
    for out in deck.outputs:
        if out.outType == "CSV":
            out.write(deck,problem)

def write_vtk(deck,problem):
    print ""

# Start the function __main__ at __init__ call
if __name__ == "__main__":
    main(sys.argv[1:])
