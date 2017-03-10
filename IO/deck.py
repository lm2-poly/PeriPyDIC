# -*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrick.diehl@polymtl.ca
import yaml
import os.path
import geometry
import sys
import util.condition
import IO.output

## Class handeling the input of the yaml file and storing the values
class PD_deck():
    
    ## Constructor
    # Reads the configuration in the yaml file and stores the values
    # @param inputFile The path to the yaml file with the configuration
    def __init__(self,inputFile):
            if not os.path.exists(inputFile):
                print "Error: Could not find " + inputFile
                sys.exit(1)
            else:
                with open(inputFile,'r') as f:
                    self.doc = yaml.load(f)
                    
                    if not "Material" in self.doc:
                        print "Error: Specify a material tag in your yaml"
                        sys.exit(1)
                    else:
                        if not "Type" in self.doc["Material"]:
                            print "Error: No type tag found"
                            sys.exit(1)
                        else:
                            ## Type of the material 
                            self.material_type = self.doc["Material"]["Type"]
                            if self.material_type == "Elastic":
                                if not "E_Modulus" in self.doc["Material"]:
                                    print "Error: No E_Modulus tag found"
                                    sys.exit(1)
                                else:
                                    ## Young's modulus of the material
                                    self.e_modulus = float(self.doc["Material"]["E_Modulus"])
                            elif self.material_type == "Viscoelastic":
                                if not "Relax_Modulus" in self.doc["Material"]:
                                    print "Error: No Relax_Modulus tag found"
                                    sys.exit(1)
                                else:
                                    ## Relaxation modulus of the material
                                    self.relax_modulus = self.doc["Material"]["Relax_Modulus"]
                                if not "Relax_Time" in self.doc["Material"]:
                                    print "Error: No Relax_Time tag found"
                                    sys.exit(1)
                                else:
                                    ## Relaxiation times
                                    self.relax_time = self.doc["Material"]["Relax_Time"]
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
                            ## Dimension of the problem
                            self.dim = int(self.doc["Discretization"]["Dim"])
                        if not "Final_Time" in self.doc["Discretization"]:
                            print "Error: No Final_Time tag found"
                            sys.exit(1)
                        else:
                            ## Final time of the problem
                            self.final_time = float(self.doc["Discretization"]["Final_Time"])
                        if not "Time_Steps" in self.doc["Discretization"]:
                            print "Error: No Time_Steps tag found"
                            sys.exit(1)
                        else:
                            ## Amount of time steps 
                            self.time_steps = int(self.doc["Discretization"]["Time_Steps"]) + 1
                            ## Time step size
                            self.delta_t = float(self.final_time  / (self.time_steps-1)) 
                        if not "Horizon_Factor_m_value" in self.doc["Discretization"]:
                            print "Error: No Horizon_Factor tag found"
                            sys.exit(1)
                        else:
                            ## The m value of the horizon factor
                            self.horizon_factor_m_value = float(self.doc["Discretization"]["Horizon_Factor_m_value"])
                        if not "Influence_Function" in self.doc["Discretization"]:
                            print "Error: Influence_Function tag found"
                            sys.exit(1)
                        else:
                            ## The influence function
                            self.influence_function = float(self.doc["Discretization"]["Influence_Function"])
                        if not ("File") in self.doc["Discretization"]:
                            print "Error: No File tag found"
                            sys.exit(1)
                        ## Object for handling the disctre nodes
                        self.geometry = geometry.Geometry()
                        self.geometry.readNodes(self.dim,self.doc["Discretization"]["File"]["Name"])
                        ## The ninimal nodal spacing in x direction
                        self.delta_x = self.geometry.getMinDist(1)
                        ## Amount of nodes in x direction
                        self.num_nodes_x = len(self.geometry.pos_x)
                        
                        if not "Boundary" in self.doc:
                            print "Error: No Boundary tag found"
                            sys.exit(1)
                        else:
                            if not "Condition" in self.doc["Boundary"]:
                                print "Error: No Condition tag found"
                                sys.exit(1)
                            else: 
                                ## List of all conditions specified in the configuration file
                                self.conditions = []
                                for i in range(0,len(self.doc["Boundary"]["Condition"]["Value"])):
                                    self.conditions.append(util.condition.ConditionFromFile(self.doc["Boundary"]["Condition"]["Type"][i],self.doc["Boundary"]["Condition"]["File"][i],self.doc["Boundary"]["Condition"]["Value"][i],self.geometry.volumes))
                            if not "Shape" in self.doc["Boundary"]:
                                print "Error: No Shape tag found"
                                sys.exit(1)                       
                            else:
                                ## Type of the shape, e.g. Ramp
                                self.shape_type = self.doc["Boundary"]["Shape"]["Type"]
                                ## List of the values for speciying the shape
                                self.shape_values = self.doc["Boundary"]["Shape"]["Values"]
                                
                        if "Output" in self.doc:
                            if  "CSV" in self.doc["Output"]:
                                if not "Type" in self.doc["Output"]["CSV"]:
                                    print "Error: No Type tag found"
                                    sys.exit(1)
                                elif not "File" in self.doc["Output"]["CSV"]:
                                    print "Error: No File tag found"
                                    sys.exit(1)
                                else:
                                    ## List of all outputs specified in the configuration file
                                    self.outputs = []
                                    for i in range(0,len(self.doc["Output"]["CSV"]["File"])):
                                        self.outputs.append(IO.output.OutputCSV("CSV",self.doc["Output"]["CSV"]["Type"][i],self.doc["Output"]["CSV"]["File"][i]))
                                      
