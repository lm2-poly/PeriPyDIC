# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 19:24:33 2015

@author: ilyass.tabiai@gmail.com
@author: rolland.delorme@gmail.com
@author: diehl@ins.uni-bonn.de
"""
import sys
import getopt
import logging
import pdb
import numpy as np
import random
from deck import PD_deck
from problem import PD_problem
from elastic_dic import elastic_material_dic
logger = logging.getLogger(__name__)

def main(argv):
    """
    Main	
    """
    types = ['elastic', 'elastic_dic']
    path = ''
    output = ''
    materialType = -1
    helpText = "__init__.py -i <inputfile> -o <outputfile> -t <type> \n" \
        "Help text"
    print sys.argv
    if len(sys.argv) != 7:
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
    if materialType not in types:
        print("Error: Only elastic and elastic_dic types are supported")
        sys.exit(1)
        
    if materialType == "elastic":
        data = PD_deck(path)
        problem = PD_problem(data)
        problem.quasi_static_solver( problem.x, data )
        problem.strain_energy_from_force( data)
        pdb.set_trace()
        
    elif materialType == "elastic_dic":
        data = PD_deck(path)
        problem = PD_problem( data )
        problem.read_csv_from_dic("data/tracking_results.csv")
        print "Experimental displacements:"
        print problem.exp_displacement
        print "Initial positions:"
        print problem.exp_init_positions
        exp_w_all = []
        for t in range(0, 13):
            elastic_dic = elastic_material_dic( data, problem, t)
            exp_w_all.append( elastic_dic.exp_W )
        problem.plot_energy(exp_w_all,problem.exp_times,problem.exp_init_positions,output)
        #problem.generate_neighborhood_matrix(data, [1, 5, 15, 25] )
        
        

#Start the function __main__ at __init__ call
if __name__ == "__main__":
    main(sys.argv[1:])





