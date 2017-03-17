import vtk
import numpy as np

class VTK_writer():
    
    def __init__(self,path,types,slice_length):
        ## Path for the output
        self.path = path
        ## Types of the attributes
        self.types = types
        ## Slice for the timesteps
        self.slice_length = slice_length
        
    def write_data(self,deck,problem):
        
        num_nodes = len(problem.y[0])
        for t in range(1,num_nodes,self.slice_length):
            writer = vtk.vtkXMLUnstructuredGridWriter()
            writer.SetFileName(self.path+"output_"+str(t)+".vtu")
            grid = vtk.vtkUnstructuredGrid()
            points = vtk.vtkPoints()
            points.SetNumberOfPoints(num_nodes)
            points.SetDataTypeToDouble()
            for i in range(num_nodes):
                if deck.dim == 1:
                    points.InsertPoint(i,problem.y[t][i],0.,0.)
                if deck.dim == 2:
                    points.InsertPoint(i,problem.y[t][i],problem.y[t][num_nodes+i],0.)
                if deck.dim == 2:
                    points.InsertPoint(i,problem.y[t][i],problem.y[t][num_nodes+i],problem.y[t][2*num_nodes+i])
                grid.SetPoints(points)
                dataOut = grid.GetPointData()    
                for out_type in self.types:
                    
                    if out_type == "Displacement":
                        array = vtk.vtkDoubleArray()
                        array.SetName("Displacement")
                        array.SetNumberOfComponents(deck.dim)
                        array.SetNumberOfTuples(num_nodes)
                        for i in range(num_nodes):
                            array.SetTuple1(i,abs(problem.y[t][i] - deck.geometry.nodes[i][0]))
                    dataOut.AddArray(array)
            
        
            writer.SetInputData(grid)
        
            writer.GetCompressor().SetCompressionLevel(0)
            writer.SetDataModeToAscii()
            writer.Write()
        
        
