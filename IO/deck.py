# -*- coding: utf-8 -*-
"""
Created on Tues Mar 07 10:00:00 2017

@author: ilyass.tabiai@polymtl.ca
@author: rolland.delorme@polymtl.ca
@author: patrick.diehl@polymtl.ca
"""

import yaml
import os.path
import geometry
import sys
import util.condition

class PD_deck():
    
    conditions = []
    material_type = ""
    e_modulus = 0.0
    relax_modulus = []
    relax_time = []
    num_nodes_x = 0.0
    delta_x = 0.0
    horizon_factor = 0.0
    time_steps = 0
    e_modulus = 0.0
    delta_t = 0.0
    shape_type = "None"
    shape_values = []
    
    def __init__(self,inputFile):
            if not os.path.exists(inputFile):
                print "Error: Could not find " + inputFile
                sys.exit(1)
            else:
                with open(inputFile,'r') as f:
                    self.doc = yaml.load(f)
                    if not "Material" in self.doc:
                        print "Error: Specific a material tag in your yaml"
                        sys.exit(1)
                    if not "Type" in self.doc["Material"]:
                        print "Error: No type tag found"
                        sys.exit(1)
                    else:
                        self.material_type = self.doc["Material"]["Type"]
                        if self.material_type == "Elastic":
                            if not "E_Modulus" in self.doc["Material"]:
                                print "Error: No E_Modulus tag found"
                                sys.exit(1)
                            else:
                                self.e_modulus = float(self.doc["Material"]["E_Modulus"])
                        elif self.material_type == "Viscoelastic":
                            if not "Relax_Modulus" in self.doc["Material"]:
                                print "Error: No Relax_Modulus tag found"
                                sys.exit(1)
                            else:
                                self.relax_modulus = self.doc["Material"]["Relax_Modulus"]
                                print self.relax_modulus
                            if not "Relax_Time" in self.doc["Material"]:
                                print "Error: No Relax_Time tag found"
                                sys.exit(1)
                            else:
                                self.relax_time = self.doc["Material"]["Relax_Time"]
                                print self.relax_time
                        else:
                            print "Error in deck.py: Material type unknown, please use Elastic or Viscoelastic"
                            sys.exit(1)
                    if not "Discretization" in self.doc:
                        print "Error: Specific a discretization tag in your yaml"
                        sys.exit(1)
                    else:
                        if not "Dim" in self.doc["Discretization"]:
                            print "Error: No Dim tag found"
                            sys.exit(1)
                        else:
                            self.dim = int(self.doc["Discretization"]["Dim"])
                        if not "Final_Time" in self.doc["Discretization"]:
                            print "Error: No Final_Time tag found"
                            sys.exit(1)
                        else:
                            self.final_time = float(self.doc["Discretization"]["Final_Time"])
                        if not "Time_Steps" in self.doc["Discretization"]:
                            print "Error: No Time_Steps tag found"
                            sys.exit(1)
                        else:
                            self.time_steps = int(self.doc["Discretization"]["Time_Steps"]) + 1
                            self.delta_t = float(self.final_time  / (self.time_steps-1)) 
                        if not "Horizon_Factor" in self.doc["Discretization"]:
                            print "Error: No Horizon_Factor tag found"
                            sys.exit(1)
                        else:
                            self.horizon_factor = float(self.doc["Discretization"]["Horizon_Factor"])
                        if not "Influence_Function" in self.doc["Discretization"]:
                            print "Error: Influence_Function tag found"
                            sys.exit(1)
                        else:
                            self.influence_function = float(self.doc["Discretization"]["Influence_Function"])
                        if not( "File") in self.doc["Discretization"]:
                            print "Error: No File tag found"
                            sys.exit(1)
                        self.geometry = geometry.Geometry()
                        self.geometry.readNodes(self.dim,self.doc["Discretization"]["File"]["Name"])
                        self.delta_x = self.geometry.getMinDist(1)
                        for i in range(0,len(self.doc["Boundary"]["Condition"]["Value"])):
                            self.conditions.append(util.condition.ConditionFromFile(self.doc["Boundary"]["Condition"]["Type"],self.doc["Boundary"]["Condition"]["File"][i],self.doc["Boundary"]["Condition"]["Value"][i],self.geometry.volumes))
                        self.num_nodes_x = len(self.geometry.pos_x)
                        if "Shape" in self.doc["Boundary"]:
                            self.shape_type = self.doc["Boundary"]["Shape"]["Type"]
                            self.shape_values = self.doc["Boundary"]["Shape"]["Values"]
                    if not "Solver" in self.doc:
                        print "Error: No solver tag found"
                        sys.exit(1)
                    self.solver_symmetry = self.doc["Solver"]["Symmetry"]
                        
                            
                                
                                
                                

    

                    
