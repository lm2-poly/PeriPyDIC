# -*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrick.diehl@polymtl.ca
import vtk

## Class for writing the attributes to an vtk unstructured grid
class vtk_writer():
    ## Constructor
    # @param path The path for the output files
    # @param types The attributes to attend to the output file
    # @param slice_length The methods writes all time steps from the first one to the last one
    #   with slice_length as the step width
    def __init__(self, path, types, slice_length):
        ## Path for the output
        self.path = path
        ## Types of the attributes
        self.types = types
        ## Slice for the timesteps
        self.slice_length = slice_length
    ## Write the data to the harddrive
    # @param deck The input deck
    # @param problem The PD or DIC problem
    def write_data(self, deck, problem):
        num_nodes = len(problem.y[0])
        for t in range(1,num_nodes, self.slice_length):
            writer = vtk.vtkXMLUnstructuredGridWriter()
            writer.SetFileName(self.path+"output_"+str(t)+".vtu")
            grid = vtk.vtkUnstructuredGrid()
            points = vtk.vtkPoints()
            points.SetNumberOfPoints(num_nodes)
            points.SetDataTypeToDouble()
            for i in range(num_nodes):
                if deck.dim == 1:
                    points.InsertPoint(i, problem.y[t][i], 0., 0.)
                if deck.dim == 2:
                    points.InsertPoint(i, problem.y[t][i], problem.y[t][num_nodes+i], 0.)
                if deck.dim == 2:
                    points.InsertPoint(i, problem.y[t][i], problem.y[t][num_nodes+i], problem.y[t][2*num_nodes+i])
                grid.SetPoints(points)
                dataOut = grid.GetPointData()
                for out_type in self.types:
                    if out_type == "Displacement":
                        array = vtk.vtkDoubleArray()
                        array.SetName("Displacement")
                        array.SetNumberOfComponents(deck.dim)
                        array.SetNumberOfTuples(num_nodes)
                        for i in range(num_nodes):
                            if deck.dim == 1:
                                array.SetTuple1(i,abs(problem.y[t][i] - deck.geometry.nodes[i][0]))
                            if deck.dim == 2:
                                dx = abs(problem.y[t][i] - deck.geometry.nodes[i][0])
                                dy = abs(problem.y[t][num_nodes + i] - deck.geometry.nodes[i][1])
                                array.SetTupl2(i, dx , dy)
                            if deck.dim == 3:
                                dx = abs(problem.y[t][i] - deck.geometry.nodes[i][0])
                                dy = abs(problem.y[t][num_nodes + i] - deck.geometry.nodes[i][1])
                                dy = abs(problem.y[t][2*num_nodes + i] - deck.geometry.nodes[i][2])
                                array.SetTupl3(i, dx, dy, dz)
                    dataOut.AddArray(array)
            writer.SetInputData(grid)
            writer.GetCompressor().SetCompressionLevel(0)
            writer.SetDataModeToAscii()
            writer.Write()
