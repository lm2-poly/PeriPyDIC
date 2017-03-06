import yaml
import os.path
import geometry
import sys
import util.condition

class Deck():
    
    conditions = []
    
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
                    if not "E_Modulus" in self.doc["Material"]:
                            print "Error: No E_Modulus tag found"
                            sys.exit(1)
                    else:
                            self.e_modulus = float(self.doc["Material"]["E_Modulus"])
                        
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
                            self.time_steps = float(self.doc["Discretization"]["Time_Steps"])
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
                        if not( "Bar") in self.doc["Discretization"] and not( "File") in self.doc["Discretization"]:
                            print "Error: No Bar or File tag found"
                            sys.exit(1)
                        if "Bar" in self.doc["Discretization"] and "File" in self.doc["Discretization"]:
                            print "Error: Bar tag and File tag found. Only one type of discretization is allowd"
                        else:
                            self.geometry = geometry.Geometry()
                            if "Bar" in self.doc["Discretization"]:
                                self.geometry.generateGrid(self.dim,self.doc["Discretization"]["Bar"],self.horizon_factor)
                                self.number_of_nodes_x = int(self.doc["Discretization"]["Bar"]["Nodes_X"])
                                self.surface_X = float(self.doc["Discretization"]["Bar"]["Surface_X"])
                                self.length_x = float(self.doc["Discretization"]["Bar"]["Length_X"])
                            else:
                                self.geometry.readNodes(self.dim,self.doc["Discretization"]["File"]["Name"])
                                self.conditions.append(util.condition.Condition(self.doc["Boundary"]["Condition"]))
                                print self.geometry.pos_x
                                print self.geometry.volumes
                                print self.conditions[0].id
                                
                                
                                

    

                    
