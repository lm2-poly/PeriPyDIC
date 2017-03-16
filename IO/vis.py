import vtk

class VTK_writer():
    
    def __init__(self,input_file):
        self.writer = vtk.vtkXMLUnstructuredGridWriter()
        self.writer.SetFileName(input_file)
        self.grid = vtk.vtkUnstructuredGrid()
        
    def write_positions(self,y):
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(numPoints)
        points.SetDataTypeToDouble()