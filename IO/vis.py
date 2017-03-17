import vtk

class VTK_writer():
    
    def __init__(self,input_file,types,slice_length):
        ## Path for the output
        self.input_file = input_file
        ## Types of the attributes
        self.types = types
        ## Slice for the timesteps
        self.slice_length = slice_length
        
    def write_data(delf,deck,problem):
        
        for t in range(0,len(y),self.slice_length):
            print t
        #self.writer = vtk.vtkXMLUnstructuredGridWriter()
        #self.writer.SetFileName(input_file)
        #self.grid = vtk.vtkUnstructuredGrid()
        
        
        
    def write_positions(self,y):
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(numPoints)
        points.SetDataTypeToDouble()
