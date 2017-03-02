# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 19:24:33 2015

@author: ilyass.tabiai@polymtl.ca
@author: rolland.delorme@polymtl.ca
@author: patrick.diehl@polymtl.ca
"""
import sys
import getopt
import logging
import pdb
import numpy as np
import random
import IO.PD_deck
from problem import PD_problem
from dic.elastic_dic import elastic_material_dic
Logger = logging.getLogger(__name__)


def main(argv):
    """
    Main
    """
    types = ['elastic', 'elastic_dic', 'viscoelastic']
    path = ''
    output = ''
    materialType = -1
    helpText = "__init__.py -i <inputfile> -o <outputfile> -t <type> \n" \
        "Help text"
    
    if len(sys.argv) != 5:
        print(helpText)
        sys.exit(1)

    try:
        opts, args = getopt.getopt(
            argv, "hi:o:t:", ["ifile=", "ofile=", "type="])
    except getopt.GetoptError:
        print(helpText)
        sys.exit(0)

    for opt, arg in opts:
        if opt == '-h':
            print(helpText)
            sys.exit(0)
        elif opt in ("-i", "--ifile"):
            path = arg
        elif opt in ("-o", "--ofile"):
            output = arg
        elif opt in ("-t", "--type"):
            materialType = arg
    
    data = IO.PD_deck.PD_deck(path)
    
    if data.Material_Flag == "ELASTIC":
#        materialType == "elastic"
        simulation(data)

    elif data.Material_Flag == "VISCOELASTIC":
#        materialType == "viscoelastic"
        simulation(data)
    
#    elif materialType == "elastic_dic":
#        data = IO.PD_deck.PD_deck(path)
#        problem = PD_problem(data)
#        problem.read_csv_from_dic("data/tracking_results.csv")
#        print "Experimental displacements:"
#        print problem.exp_displacement
#        print "Initial positions:"
#        print problem.exp_init_positions
#        exp_w_all = []
#        for t_n in range(0, data.N_Steps_t):
#            elastic_dic = elastic_material_dic(data, problem, t_n)
#            exp_w_all.append(elastic_dic.exp_W)
#        problem.plot_energy(
#            exp_w_all,
#            problem.exp_times,
#            problem.exp_init_positions,
#            output + "/elastic_dic_")
    else:
        print("Error: Material type not supported")
        sys.exit(1)

def simulation(data):
    problem = PD_problem(data)
    problem.quasi_static_solver(problem.x, data)
    problem.strain_center_bar( data )
    print np.around(problem.strain,decimals=6)
#   problem.strain_energy_from_force(data)
#   problem.plot_energy(
#       problem.remove_additional_points(data, problem.strain_energy_from_force),
#       data.time_steps,
#       problem.x[:len(problem.x)-2],
#       output + "/viscoelastic_")

# Start the function __main__ at __init__ call
if __name__ == "__main__":
    main(sys.argv[1:])
