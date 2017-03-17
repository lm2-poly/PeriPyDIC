#-*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrick.diehl@polymtl.ca
import sys
import getopt
import IO.deck
import problem
import numpy as np
import pkgutil

vtk_loader = pkgutil.find_loader('vtk')
found_vtk = vtk_loader is not None 
if found_vtk == True:
    import IO.vis
else:
    print "VTK output not supported. Please install PyVTK for that"
    
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
    solver = problem.PD_problem(deck)
    
    if deck.dim == 1:
        initialVector = deck.geometry.nodes
    if deck.dim == 2:
        initialVector = np.concatenate((deck.geometry.nodes[:,0], deck.geometry.nodes[:,1]), axis=0)
    if deck.dim == 3:
        initialVector = np.concatenate((deck.geometry.nodes[:,0] ,np.concatenate((deck.geometry.nodes[:,1], deck.geometry.nodes[:,2]), axis=0)),axis=0)
  
    x_0 = solver.random_initial_guess( initialVector, deck )
    solver.quasi_static_solver(x_0, deck)
    solver.strain_center_bar( deck )
    writeCSV(deck,solver)
    if found_vtk == True and deck.vtk_output == True:
       deck.vtk_writer.write_data(deck,solver)
    print "delta_x =" , deck.delta_x
    print "Horizon =" , solver.neighbors.horizon
    print "Strain = " , np.around(solver.strain,decimals=6)
    #print "Nodes positions = "
    #print solver.y

def writeCSV(deck,problem):
    for out in deck.outputs:
        if out.outType == "CSV":
            out.write(deck,problem)

def write_vtk(deck,problem):
    print ""
    
# Start the function __main__ at __init__ call
if __name__ == "__main__":
    main(sys.argv[1:])
