# -*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrickdiehl@lsu.edu
import yaml
import os.path
from . import geometry
import sys
from . import output
from ..util import condition
from . import vis
from . import dic
import numpy as np

## Class handeling the input of the yaml file and storing the values
class PD_deck():

    ## Constructor
    # Reads the configuration in the yaml file and stores the values
    # @param inputFile The path to the yaml file with the configuration
    def __init__(self,inputFile):
            if not os.path.exists(inputFile):
                print ("Error: Could not find " + inputFile)
                sys.exit(1)
            else:
                with open(inputFile,'r') as f:
                    ## Container of the tags parsed from the yaml file
                    self.doc = yaml.load(f, Loader=yaml.FullLoader)
                    ## Safety factor for the computation of the radius
                    self.safety_factor = 1.001
                    ## Number of threads
                    self.num_threads = 1

                    if not "Discretization" in self.doc:
                        print ("Error: Specific a Discretization tag in your yaml")
                        sys.exit(1)
                    else:
                        if not "Dim" in self.doc["Discretization"]:
                            print ("Error: No Dim tag found")
                            sys.exit(1)
                        else:
                            ## Dimension of the problem
                            self.dim = int(self.doc["Discretization"]["Dim"])
                            if self.dim == 2:
                                if not "Type" in self.doc["Discretization"]:
                                    print ("Error: No Type tag found for 2D problem. Choose Plane_Stress or Plane_Strain")
                                    sys.exit(1)
                                else:
                                    ## Type of 2D problem
                                    self.type2d = self.doc["Discretization"]["Type"]
                        if not "Final_Time" in self.doc["Discretization"]:
                            print ("Error: No Final_Time tag found")
                            sys.exit(1)
                        else:
                            ## Final time of the problem
                            self.final_time = float(self.doc["Discretization"]["Final_Time"])
                        if not "Time_Steps" in self.doc["Discretization"]:
                            print ("Error: No Time_Steps tag found")
                            sys.exit(1)
                        else:
                            ## Amount of time steps
                            self.time_steps = int(self.doc["Discretization"]["Time_Steps"]) + 1
                            ## Time step size
                            self.delta_t = float(self.final_time  / (self.time_steps-1))
                        if not "Horizon_Factor_m_value" in self.doc["Discretization"]:
                            print ("Error: No Horizon_Factor tag found")
                            sys.exit(1)
                        else:
                            ## "m" value of the horizon factor
                            self.horizon_factor_m_value = float(self.doc["Discretization"]["Horizon_Factor_m_value"])
                        if not "Influence_Function" in self.doc["Discretization"]:
                            print ("Error: Influence_Function tag found")
                            sys.exit(1)
                        else:
                            ## Influence function
                            self.influence_function = self.doc["Discretization"]["Influence_Function"]
                        if "Safety_Factor" in self.doc["Discretization"]:
                            self.safety_factor = float(self.doc["Discretization"]["Safety_Factor"])
                        if not ("File") in self.doc["Discretization"]:
                            print ("Error: No File tag found")
                            sys.exit(1)
                        ## Object for handling the discrete nodes
                        self.geometry = geometry.Geometry()
                        self.geometry.readNodes(self.dim,self.doc["Discretization"]["File"]["Name"])
                        ## The minimal nodal spacing
                        self.delta_X = self.geometry.getMinDist()
                        ## Amount of nodes
                        self.num_nodes = self.geometry.amount
                        if "Boundary" in self.doc:
                            if not "Condition" in self.doc["Boundary"]:
                                print ("Error: No Condition tag found")
                                sys.exit(1)
                            else:
                                ## List of all conditions specified in the configuration file
                                self.conditions = []
                                for i in range(0,len(self.doc["Boundary"]["Condition"]["Value"])):
                                    self.conditions.append(condition.ConditionFromFile(self.doc["Boundary"]["Condition"]["Type"][i],self.doc["Boundary"]["Condition"]["File"][i],self.doc["Boundary"]["Condition"]["Value"][i],self.geometry.volumes,self.doc["Boundary"]["Condition"]["Direction"][i],self.doc["Boundary"]["Condition"]["Shape"][i]))
                            if not "Shape" in self.doc["Boundary"]:
                                print ("Error: No Shape tag found")
                                sys.exit(1)
                            else:
                                ## Type of the shape, e.g. Ramp
                                self.shape_type = self.doc["Boundary"]["Shape"]["Type"]
                                ## List of the values for specifying the shape
                                self.shape_values = self.doc["Boundary"]["Shape"]["Values"]

                    if not "Material" in self.doc:
                        print ("Error: Specify a material tag in your yaml")
                        sys.exit(1)
                    else:
                        if not "Type" in self.doc["Material"]:
                            print ("Error: No Type tag found")
                            sys.exit(1)
                        else:
                            ## Type of the material
                            self.material_type = self.doc["Material"]["Type"]
                            if self.material_type == "Elastic":
                                if self.dim == 1:
                                    if not "Young_Modulus" in self.doc["Material"]:
                                        print ("Error: No Young_Modulus tag found")
                                        sys.exit(1)
                                    else:
                                        ## Young modulus of the material
                                        self.young_modulus = float(self.doc["Material"]["Young_Modulus"])
                                if self.dim >= 2:
                                    if not "Bulk_Modulus" in self.doc["Material"]:
                                        print ("Error: No Bulk_Modulus tag found")
                                        sys.exit(1)
                                    else:
                                        ## Bulk modulus of the material
                                        self.bulk_modulus = float(self.doc["Material"]["Bulk_Modulus"])
                                    if not "Shear_Modulus" in self.doc["Material"]:
                                        print ("Error: No Shear_Modulus tag found")
                                        sys.exit(1)
                                    else:
                                        ## Shear modulus of the material
                                        self.shear_modulus = float(self.doc["Material"]["Shear_Modulus"])
                            elif self.material_type == "Viscoelastic":
                                if self.dim == 1:
                                    if not "Relax_Modulus" in self.doc["Material"]:
                                        print ("Error: No Relax_Modulus tag found")
                                        sys.exit(1)
                                    else:
                                        ## Relaxation modulus of the material
                                        self.relax_modulus = self.doc["Material"]["Relax_Modulus"]
                                if self.dim >= 2:
                                    if not "Relax_Bulk_Modulus" in self.doc["Material"]:
                                        print ("Error: No Relax_Bulk_Modulus tag found")
                                        sys.exit(1)
                                    else:
                                        ## Relaxation Bulk modulus of the material
                                        self.relax_bulk_modulus = np.asarray(self.doc["Material"]["Relax_Bulk_Modulus"], dtype=np.float64)
                                    if not "Relax_Shear_Modulus" in self.doc["Material"]:
                                        print ("Error: No Relax_Shear_Modulus tag found")
                                        sys.exit(1)
                                    else:
                                        ## Relaxation Shear modulus of the material
                                        self.relax_shear_modulus = np.asarray(self.doc["Material"]["Relax_Shear_Modulus"], dtype=np.float64)
                                if not "Relax_Time" in self.doc["Material"]:
                                    print ("Error: No Relax_Time tag found")
                                    sys.exit(1)
                                else:
                                    ## Relaxation times
                                    self.relax_time = np.asarray(self.doc["Material"]["Relax_Time"], dtype=np.float64)
                            else:
                                print ("Error in deck.py: Material type unknown, please use Elastic or Viscoelastic")
                                sys.exit(1)
                        ## List of all outputs specified in the configuration file
                        self.outputs = []
                        if "Output" in self.doc:
                            if  "CSV" in self.doc["Output"]:
                                if not "Type" in self.doc["Output"]["CSV"]:
                                    print ("Error: No Type tag found")
                                    sys.exit(1)
                                elif not "File" in self.doc["Output"]["CSV"]:
                                    print ("Error: No File tag found")
                                    sys.exit(1)
                                else:
                                    for i in range(0,len(self.doc["Output"]["CSV"]["File"])):
                                        self.outputs.append(output.OutputCSV("CSV",self.doc["Output"]["CSV"]["Type"][i],self.doc["Output"]["CSV"]["File"][i]))
                            if "VTK" in self.doc["Output"]:
                                if not "Path" in self.doc["Output"]["VTK"]:
                                    print ("Error: No Path tag found in VTK")
                                    sys.exit(1)
                                elif not "Type" in self.doc["Output"]["VTK"]:
                                    print ("Error: No Type tag found in VTK")
                                    sys.exit(1)
                                elif not "Slice" in self.doc["Output"]["VTK"]:
                                    print ("Error: No Slice tag found in VTK")
                                    sys.exit(1)
                                else:
                                    ## Visualization ToolKit (VTK) writer
                                    self.vtk_writer = vis.vtk_writer(self.doc["Output"]["VTK"]["Path"],self.doc["Output"]["VTK"]["Type"],self.doc["Output"]["VTK"]["Slice"])
                                    if self.vtk_writer.vtk_enabled == False:
                                        print ("Warning: VTK found, but no PyVTK is found, so there will be no output written.")
                            else:
                                self.vtk_writer = vis.vtk_writer()
                        else:
                            self.vtk_writer = vis.vtk_writer()

                        if not "Solver" in  self.doc:
                            print ("Error: No Solver tag found")
                            sys.exit(1)
                        else:
                            if not "Max_Iteration" in self.doc["Solver"]:
                                print ("Error: No Max_Iteration tag in Solver found")
                                sys.exit(1)
                            else:
                                ## Maximum number iteration
                                self.solver_max_it = self.doc["Solver"]["Max_Iteration"]
                            if not "Tolerance" in self.doc["Solver"]:
                                print ("Error: No Tolerance tag in Solver found")
                                sys.exit(1)
                            else:
                                ## Absolute tolerance of the solver
                                self.solver_tolerance = float(self.doc["Solver"]["Tolerance"])
                            if not "Jacobian_Perturbation" in self.doc["Solver"]:
                                print ("Error: No Jacobian_Perturbation tag in Solver found")
                                sys.exit(1)
                            else:
                                ## Perturbation factor for the Jacobian matrix
                                self.solver_perturbation = float(self.doc["Solver"]["Jacobian_Perturbation"])

                        if "Parallel" in self.doc:
                            if "Threads" in self.doc["Parallel"]:
                                self.num_threads = int(self.doc["Parallel"]["Threads"])

class DIC_deck():

    ## Constructor
    # Reads the configuration in the yaml file and stores the values
    # @param inputFile The path to the yaml file with the configuration
    def __init__(self,inputFile):
            if not os.path.exists(inputFile):
                print ("Error: Could not find " + inputFile)
                sys.exit(1)
            else:
                with open(inputFile,'r') as f:
                    ## Container of the tags parsed from the yaml file
                    self.doc = yaml.load(f, Loader=yaml.FullLoader)
                    if not "Material" in self.doc:
                        print ("Error: Specify a Material tag in your yaml")
                        sys.exit(1)
                    else:
                        if not "Type" in self.doc["Material"]:
                            print ("Error: No type tag found")
                            sys.exit(1)
                        else:
                            ## Type of the material
                            self.material_type = self.doc["Material"]["Type"]
                            if self.material_type == "Elastic":
                                if "Young_Modulus" in self.doc["Material"]:
                                    ## Young modulus of the material
                                    self.young_modulus = float(self.doc["Material"]["Young_Modulus"])
                                if "Bulk_Modulus" in self.doc["Material"]:
                                    ## Bulk modulus of the material
                                    self.bulk_modulus = float(self.doc["Material"]["Bulk_Modulus"])
                                if "Shear_Modulus" in self.doc["Material"]:
                                    ## Shear modulus of the material
                                    self.shear_modulus = float(self.doc["Material"]["Shear_Modulus"])

                            elif self.material_type == "Viscoelastic":
                                if not "Relax_Modulus" in self.doc["Material"]:
                                    print ("Error: No Relax_Modulus tag found")
                                    sys.exit(1)
                                else:
                                    ## Relaxation modulus of the material
                                    self.relax_modulus = self.doc["Material"]["Relax_Modulus"]
                                if not "Relax_Time" in self.doc["Material"]:
                                    print ("Error: No Relax_Time tag found")
                                    sys.exit(1)
                                else:
                                    ## Relaxation times
                                    self.relax_time = self.doc["Material"]["Relax_Time"]
                            else:
                                print ("Error in deck.py: Material type unknown, please use Elastic or Viscoelastic")
                                sys.exit(1)

                        if not "Discretization" in self.doc:
                            print ("Error: Specify a Discretization tag in your yaml")
                            sys.exit(1)
                        else:
                            ## Safety factor for the computation of the radius
                            self.safety_factor = 1.001
                            if not "Horizon_Factor_m_value" in self.doc["Discretization"]:
                                print ("Error: No Horizon_Factor_m_value tag found")
                                sys.exit(1)
                            else:
                                ## "m" value of the horizon factor
                                self.horizon_factor_m_value = self.doc["Discretization"]["Horizon_Factor_m_value"]

                            if not "Influence_Function" in self.doc["Discretization"]:
                                print ("Error: No Influence_Function tag found")
                                sys.exit(1)
                            else:
                                ## Influence function
                                self.influence_function = self.doc["Discretization"]["Influence_Function"]
                            if "Saftety_Factor" in self.doc["Discretization"]:
                                self.safety_factor = float(self.doc["Discretization"]["Safety_Factor"])

                        if not "Data" in self.doc:
                            print ("Error: Specify a Data tag in your yaml")
                            sys.exit(1)
                        else:
                            if not "Dimension" in self.doc["Data"]:
                                print ("Error: No Dimension tag found")
                                sys.exit(1)
                            else:
                                ## The dimension of the input data
                                self.dim = self.doc["Data"]["Dimension"]
                                if self.dim == 2:
                                    if not "Type" in self.doc["Data"]:
                                        print ("Error: No Type tag found for 2D problem. Choose Plane_Stress or Plane_Strain")
                                        sys.exit(1)
                                    else:
                                        ## Type of 2D problem
                                        self.type2d = self.doc["Data"]["Type"]

                            if not "Sigma" in self.doc["Data"]:
                                print ("Error: No Sigma (column of confidence in CSV file) tag found")
                                sys.exit(1)
                            else:
                                ## The column number of the confidence in VIC3D CSV file
                                self.sigma_column = self.doc["Data"]["Sigma"]

                            if not "File" in self.doc["Data"]:
                                print ("Error: Specify a File tag in your yaml")
                                sys.exit(1)
                            else:
                                if not "Name" in self.doc["Data"]["File"]:
                                    print ("Error: No Name tag found")
                                    sys.exit(1)
                                else:
                                    ## Filename of the input file
                                    self.filename = self.doc["Data"]["File"]["Name"]

                                if not "Type" in self.doc["Data"]["File"]:
                                    print ("Error: No Type tag found")
                                    sys.exit(1)
                                else:
                                    self.filetype = self.doc["Data"]["File"]["Type"]

                                if not "Path" in self.doc["Data"]["File"]:
                                    print ("Error: No Path tag found")
                                    sys.exit(1)
                                else:
                                    ## File path
                                    self.filepath = self.doc["Data"]["File"]["Path"]
                                    ## Nodes uploaded from DIC input file
                                    self.geometry = dic.DICreader2D(self)
                                    ## Amount of nodes
                                    self.num_nodes = len(self.geometry.nodes)
                                    ## Minimal nodal spacing
                                    self.delta_X = self.geometry.delta_x

                                if not "Volume" in self.doc["Discretization"] and self.filetype == "vic3d":
                                    print ("Error: No Volume tag found in Discretization which is needed VIC3D")
                                    sys.exit(1)
                                else:
                                    ## Volume of a DIC node
                                    self.dic_volume = float(self.doc["Discretization"]["Volume"])
                                

                        ## Amount of time steps
                        self.time_steps = 2
                        ## Number of threads
                        self.num_threads = 1
                        if "Output" in self.doc:
                            if "VTK" in self.doc["Output"]:
                                if not "Path" in self.doc["Output"]["VTK"]:
                                    print ("Error: No Path tag found in VTK")
                                    sys.exit(1)
                                elif not "Type" in self.doc["Output"]["VTK"]:
                                    print ("Error: No Type tag found in VTK")
                                    sys.exit(1)
                                else:
                                    ## Visualization ToolKit (VTK) writer
                                    self.vtk_writer = vis.vtk_writer(self.doc["Output"]["VTK"]["Path"],self.doc["Output"]["VTK"]["Type"],1)
                                    if self.vtk_writer.vtk_enabled == False:
                                        print ("Warning: VTK found, but no PyVTK is found, so there will be no output written.")
                        if "Parallel" in self.doc:
                            if "Threads" in self.doc["Parallel"]:
                                self.num_threads = int(self.doc["Parallel"]["Threads"])
                                
                        if "Energy" in self.doc:
                            if not "Measured Energy" in self.doc["Energy"]:
                                print ("Error: Measured Energy tag found in Energy")
                                sys.exit(1)
                            else:
                                ## Measured energy from DIC
                                self.measured_energy = float(self.doc["Energy"]["Measured Energy"])
                            if not "Nodes" in self.doc["Energy"]:
                                print ("Error: No Nodes tag found in Energy")
                                sys.exit(1)
                            else:
                                ## Nodes to comapre the energy with
                                self.nodes_compare = np.asarray(self.doc["Energy"]["Nodes"], dtype=np.int)
                                self.compare_length = len(self.nodes_compare)
                                
                        if "Solver" in  self.doc:
                            if not "Max_Iteration" in self.doc["Solver"]:
                                print ("Error: No Max_Iteration tag in Solver found")
                                sys.exit(1)
                            else:
                                ## Maximum number iteration
                                self.solver_max_it = self.doc["Solver"]["Max_Iteration"]
                            if not "Tolerance" in self.doc["Solver"]:
                                print ("Error: No Tolerance tag in Solver found")
                                sys.exit(1)
                            else:
                                ## Absolute tolerance of the solver
                                self.solver_tolerance = float(self.doc["Solver"]["Tolerance"])
                            if not "Jacobian_Perturbation" in self.doc["Solver"]:
                                print ("Error: No Jacobian_Perturbation tag in Solver found")
                                sys.exit(1)
                            else:
                                ## Perturbation factor for the Jacobian matrix
                                self.solver_perturbation = float(self.doc["Solver"]["Jacobian_Perturbation"])
