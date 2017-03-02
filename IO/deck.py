
import yaml
import os.path
import geometry

class Deck():
    
    def __init__(self,inputFile):
            if not os.path.exists(inputFile):
                print "Error: Could not find " + inputFile
            else:
                with open(inputFile,'r') as f:
                    self.doc = yaml.load(f)
                    
                    if not "Material" in self.doc:
                        print "Error: Specific a material tag in your yaml"
                        if not "Type" in self.doc["Material"]:
                            print "Error: No type tag found"
                        else:
                            self.material_type = self.doc["Material"]["Type"]
                    if not "E_Modulus" in self.doc["Material"]:
                            print "Error: No E_Modulus tag found"
                    else:
                            self.e_modulus = float(self.doc["Material"]["E_Modulus"])
                        
                    if not "Discretization" in self.doc:
                        print "Error: Specific a discretization tag in your yaml"
                    else:
                        if not "Dim" in self.doc["Discretization"]:
                            print "Error: No Dim tag found"
                        else:
                            self.dim = int(self.doc["Discretization"]["Dim"])
                        if not "Final_Time" in self.doc["Discretization"]:
                            print "Error: No Final_Time tag found"
                        else:
                            self.final_time = float(self.doc["Discretization"]["Final_Time"])
                        if not "Time_Steps" in self.doc["Discretization"]:
                            print "Error: No Time_Steps tag found"
                        else:
                            self.time_steps = float(self.doc["Discretization"]["Time_Steps"])
                        if not "Horizon_Factor" in self.doc["Discretization"]:
                            print "Error: No Horizon_Factor tag found"
                        else:
                            self.horizon_factor = float(self.doc["Discretization"]["Horizon_Factor"])
                        if not "Influence_Function" in self.doc["Discretization"]:
                            print "Error: Influence_Function tag found"
                        else:
                            self.influence_function = float(self.doc["Discretization"]["Influence_Function"])
                        if not( "Bar" or "File") in self.doc["Discretization"]:
                            print "Error: No Bar or File tag found"
                        else:
                            self.geometry = geometry.Geometry()
                            if "Bar" in self.doc:
                                self.geometry.generateGrid(self.dim,doc["Bar"],self.horizon_factor)
                                self.number_of_nodes_x = int(self.doc["Discretization"]["Nodes_X"])
                                self.surface_X = float(self.doc["Discretization"]["Surface_X"])
                                self.length_x = float(self.doc["Discretization"]["Length_X"])
                            elif "File" in self.doc:
                                self.geometry.readNodes(self.dim,doc["File"])
                                
                                
                                

    

                    
