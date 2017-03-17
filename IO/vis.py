import vtk

class VTK_writer():
    
    def __init__(self,path,types,slice_length):
        ## Path for the output
        self.path = path
        ## Types of the attributes
        self.types = types
        ## Slice for the timesteps
        self.slice_length = slice_length
        
    def write_data(delf,deck,problem):
        
        num_nodes = len(deck.y[0])
        for t in range(0,num_nodes,self.slice_length):
            writer = vtk.vtkXMLUnstructuredGridWriter()
            writer.SetFileName(input_file+"output_"+str(t)+"*.vtu")
            grid = vtk.vtkUnstructuredGrid()
        
        
        
    def write_positions(self,y):
        points = vtk.vtkPoints()
        
        points.SetDataTypeToDouble()
