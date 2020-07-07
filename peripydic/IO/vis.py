# -*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrickdiehl@lsu.edu

import pkgutil
vtk_loader = pkgutil.find_loader('vtk')
found_vtk = vtk_loader is not None
if found_vtk == True:
    import vtk
## Handles the output of the simulation results to vtk unstrucutred grids
class vtk_writer():

    if found_vtk == True:

        ## Constructor
        #
        # @param path The path where the output is written
        # @param types The simulation attributes considered as output
        # @param slice_length Every n-th time step is written
        def __init__(self,path="",types="",slice_length=-1):
            ## IS vtk enabled
            self.vtk_enabled = True
            ## Path for the output
            self.path = path
            ## Types of the attributes
            self.types = types
            ## Slice for the time steps
            self.slice_length = slice_length

        ## Writes the data to a vtk unstrucutred grid files
        #
        # @param deck The deck with the input from the yaml file
        # @param problem The problem contains the simulation results
        # @param ccm_class The results from the computation from ccm
        def write_data(self,deck,problem,ccm_class):
            num_nodes = deck.num_nodes
            for t in range(0,deck.time_steps,self.slice_length):
                writer = vtk.vtkXMLUnstructuredGridWriter()
                writer.SetFileName(self.path+"output_"+str(t)+".vtu")
                grid = vtk.vtkUnstructuredGrid()
                points = vtk.vtkPoints()
                points.SetNumberOfPoints(num_nodes)
                points.SetDataTypeToDouble()
                for i in range(0,num_nodes):
                    act = problem.y
                    if deck.dim == 1:
                        points.InsertPoint(i,act[i][0][t],0.,0.)
                    if deck.dim == 2:
                        points.InsertPoint(i,act[i][0][t],act[i][1][t],0.)
                    if deck.dim == 3:
                        points.InsertPoint(i,act[i][0][t],act[i][1][t],act[i][2][t])
                    grid.SetPoints(points)

                dataOut = grid.GetPointData()

                for out_type in self.types:

                    if out_type == "Displacement":
                        array = vtk.vtkDoubleArray()
                        array.SetName("Displacement")
                        array.SetNumberOfComponents(deck.dim)
                        array.SetNumberOfTuples(num_nodes)

                        act = problem.y

                        for i in range(num_nodes):
                            if deck.dim == 1:
                                array.SetTuple1(i,act[i][0][t] - deck.geometry.nodes[i][0])
                            if deck.dim == 2:
                                array.SetTuple2(i,act[i][0][t] - deck.geometry.nodes[i][0],act[i][1][t] - deck.geometry.nodes[i][1])
                                array.SetComponentName(0,"d_x")
                                array.SetComponentName(1,"d_y")
                            dataOut.AddArray(array)

                    if out_type == "Neighbors":
                        array = vtk.vtkIntArray()
                        array.SetName("Neighbors")
                        array.SetNumberOfComponents(1)
                        array.SetNumberOfTuples(num_nodes)

                        for i in range(num_nodes):
                            array.SetTuple1(i,len(problem.neighbors.get_index_x_family(i)))
                        dataOut.AddArray(array)

                    if out_type == "Force":
                        array = vtk.vtkDoubleArray()
                        array.SetName("Volume_Force")
                        array.SetNumberOfComponents(deck.dim)
                        array.SetNumberOfTuples(num_nodes)

                        force = problem.force_int
                        #print force
                        for i in range(num_nodes):
                            if deck.dim == 1:
                                array.SetTuple1(i,force[i][0][t])
                            if deck.dim == 2:
                                array.SetTuple2(i,force[i][0][t], force[i][1][t])
                                array.SetComponentName(0,"f_x")
                                array.SetComponentName(1,"f_y")
                            dataOut.AddArray(array)

                    if out_type == "Conditions":

                        for con in deck.conditions:
                            array = vtk.vtkIntArray()
                            array.SetName("Condition_"+con.type+"_"+str(con.value)+"_"+str(con.direction))
                            array.SetNumberOfComponents(1)
                            array.SetNumberOfTuples(num_nodes)

                            for i in range(num_nodes):
                                if i not in con.id:
                                    array.SetTuple1(i,0)
                                else:
                                    array.SetTuple1(i,1)
                        dataOut.AddArray(array)

                    if out_type == "Volume_Force":

                        force = problem.force_int
                        for con in deck.conditions:
                            if con.type == "Force":
                                result_x = 0.
                                result_y = 0.
                                result_z = 0.
                                for i in con.id:
                                    index = int(i)

                                    if deck.dim >=1:
                                        result_x += force[index][0][t] * deck.geometry.volumes[index]
                                    if deck.dim >= 2:
                                        result_y += force[index][1][t] * deck.geometry.volumes[index]
                                        array.SetComponentName(0,"f_x")
                                        array.SetComponentName(1,"f_y")
                                    if deck.dim >= 3:
                                        result_z += force[index][2][t] * deck.geometry.volumes[index]
                                        array.SetComponentName(2,"f_z")

                                array = vtk.vtkDoubleArray()
                                array.SetName("Volume_"+con.type+"_"+str(con.value)+"_"+str(con.direction))
                                array.SetNumberOfComponents(deck.dim)
                                array.SetNumberOfTuples(num_nodes)


                                for i in range(num_nodes):
                                    if i in con.id:
                                        if deck.dim ==1:
                                            array.SetTuple1(i,result_x)
                                        if deck.dim == 2:
                                            array.SetTuple2(i,result_x,result_y)
                                        if deck.dim == 3:
                                            array.SetTuple3(i,result_x,result_y,result_z)
                                    else:
                                        if deck.dim ==1:
                                            array.SetTuple1(i,0.)
                                        if deck.dim == 2:
                                            array.SetTuple2(i,0.,0.)
                                        if deck.dim == 3:
                                            array.SetTuple3(i,0.,0.,0.)
                            dataOut.AddArray(array)

                    if out_type == "Strain":
                        array = vtk.vtkDoubleArray()
                        array.SetName("Strain")
                        if deck.dim == 1:
                            array.SetNumberOfComponents(1)
                        if deck.dim == 2:
                            array.SetNumberOfComponents(3)
                        if deck.dim == 3:
                            array.SetNumberOfComponents(6)
                        array.SetNumberOfTuples(num_nodes)

                        strain = ccm_class.global_strain[:,:,1]

                        for i in range(num_nodes):
                            if deck.dim ==1:
                                array.SetComponentName(0,"epsil_xx")
                                array.SetTuple1(i,strain[i,0])
                            if deck.dim == 2:
                                xx = strain[i*deck.dim,0]
                                xy = strain[i*deck.dim,1]
                                yy = strain[i*deck.dim+1,1]
                                array.SetTuple3(i,xx,yy,xy)
                                array.SetComponentName(0,"epsil_xx")
                                array.SetComponentName(1,"epsil_yy")
                                array.SetComponentName(2,"epsil_xy")
                            if deck.dim == 3:
                                xx = strain[i*deck.dim,0]
                                xy = strain[i*deck.dim,1]
                                yy = strain[i*deck.dim+1,1]
                                yz = strain[i*deck.dim+1,2]
                                xz = strain[i*deck.dim,2]
                                zz = strain[i*deck.dim+2,2]
                                array.SetTuple6(i,xx,yy,zz,yz,xz,xy)
                                array.SetComponentName(0,"epsil_xx")
                                array.SetComponentName(1,"epsil_yy")
                                array.SetComponentName(2,"epsil_zz")
                                array.SetComponentName(3,"epsil_yz")
                                array.SetComponentName(4,"epsil_xz")
                                array.SetComponentName(5,"epsil_xy")
                        dataOut.AddArray(array)

                    if out_type == "Stress":
                        array = vtk.vtkDoubleArray()
                        array.SetName("Stress")
                        if deck.dim == 1:
                            array.SetNumberOfComponents(1)
                        if deck.dim == 2:
                            array.SetNumberOfComponents(3)
                        if deck.dim == 3:
                            array.SetNumberOfComponents(6)
                        array.SetNumberOfTuples(num_nodes)

                        stress = ccm_class.global_stress[:,:,1]

                        for i in range(num_nodes):
                            if deck.dim ==1:
                                array.SetComponentName(0,"sigma_xx")
                                array.SetTuple1(i,strain[i,0])
                            if deck.dim == 2:
                                xx = stress[i*deck.dim,0]
                                xy = stress[i*deck.dim,1]
                                yy = stress[i*deck.dim+1,1]
                                array.SetTuple3(i,xx,yy,xy)
                                array.SetComponentName(0,"sigma_xx")
                                array.SetComponentName(1,"sigma_yy")
                                array.SetComponentName(2,"sigma_xy")
                            if deck.dim == 3:
                                xx = stress[i*deck.dim,0]
                                xy = stress[i*deck.dim,1]
                                yy = stress[i*deck.dim+1,1]
                                yz = stress[i*deck.dim+1,2]
                                xz = stress[i*deck.dim,2]
                                zz = stress[i*deck.dim+2,2]
                                array.SetTuple6(i,xx,yy,zz,yz,xz,xy)
                                array.SetComponentName(0,"sigma_xx")
                                array.SetComponentName(1,"sigma_yy")
                                array.SetComponentName(2,"sigma_zz")
                                array.SetComponentName(3,"sigma_yz")
                                array.SetComponentName(4,"sigma_xz")
                                array.SetComponentName(5,"sigma_xy")
                        dataOut.AddArray(array)

                    if out_type == "Strain_DIC":
                        array = vtk.vtkDoubleArray()
                        array.SetName("Strain_DIC")
                        array.SetNumberOfComponents(3)
                        array.SetNumberOfTuples(num_nodes)

                        for i in range(num_nodes):
                            xx = deck.geometry.strain[i][0]
                            xy = deck.geometry.strain[i][2]
                            yy = deck.geometry.strain[i][1]
                            array.SetTuple3(i,xx,yy,xy)
                            array.SetComponentName(0,"epsil_xx")
                            array.SetComponentName(1,"epsil_yy")
                            array.SetComponentName(2,"epsil_xy")
                        dataOut.AddArray(array)

                    if out_type == "Strain_Error":
                        array = vtk.vtkDoubleArray()
                        array.SetName("Strain_Error")
                        array.SetNumberOfComponents(3)
                        array.SetNumberOfTuples(num_nodes)

                        strain = ccm_class.global_strain[:,:,1]

                        for i in range(num_nodes):
                            xx = abs(deck.geometry.strain[i][0] - strain[i*deck.dim,0])
                            xy = abs(deck.geometry.strain[i][2] - strain[i*deck.dim,1])
                            yy = abs(deck.geometry.strain[i][1] - strain[i*deck.dim+1,1])
                            array.SetTuple3(i,xx,yy,xy)
                            array.SetComponentName(0,"error_xx")
                            array.SetComponentName(1,"error_yy")
                            array.SetComponentName(2,"error_xy")
                        dataOut.AddArray(array)
                        
                    if out_type == "Strain_Energy":
                        array = vtk.vtkDoubleArray()
                        array.SetName("Strain_Energy")
                        array.SetNumberOfComponents(1)
                        array.SetNumberOfTuples(num_nodes)
                        
                        strain_energy = problem.strain_energy
                        
                        for i in range(num_nodes):
                            array.SetTuple1(i,strain_energy[i])
                        
                        dataOut.AddArray(array)   

                    if out_type == "Volume":
                        array = vtk.vtkDoubleArray()
                        array.SetName("Volume")
                        array.SetNumberOfComponents(1)
                        array.SetNumberOfTuples(num_nodes)

                        for i in range(num_nodes):
                            array.SetTuple1(i,deck.geometry.volumes[i])
                        
                        dataOut.AddArray(array) 

                writer.SetInputData(grid)

                writer.GetCompressor().SetCompressionLevel(0)
                writer.SetDataModeToAscii()
                writer.Write()

    else:

        ## Constructor
        #
        # @param path The path where the output is written
        # @param types The simulation attributes considered as output
        # @param slice_length Every n-th time step is written
        def __init__(self,path="",types="",slice_length=-1):
            self.vtk_enabled = False
