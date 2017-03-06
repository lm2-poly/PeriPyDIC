# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 19:24:33 2015

@author: ilyass.tabiai@gmail.com
@author: rolland.delorme@gmail.com
@author: diehl@ins.uni-bonn.de
"""

import sys
import getopt
import IO.deck

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
        deck = IO.deck.Deck(inputFile)
        
# Start the function __main__ at __init__ call
if __name__ == "__main__":
    main(sys.argv[1:])
