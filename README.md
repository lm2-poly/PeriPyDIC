# The problem

   Consider a 1D bar. The length of the bar and the number of sections in which the bar will be divided (the Peridynamic (PD) nodes are at the edges of each section) are provided by the user.
   
   A ramp loading is applied at the edges of the bar. The user provides the number of seconds after which the ramp reaches a maximal force value and remains there.
   Additional blocks are added at the edges of the block in order to apply the load. The number of blocks added at each edge of the bar is provided as the Horizon Factor.
   The user also provides the total duration of the simulation and the total number of time steps desired.
   
   The problem is then solved using an initial guess vector and the *Newton-Krylov* algorithm provided in the `scipy.optimize` package.  

# Getting started

###Dependencies

   The following python packages are **required**:
   * `numpy`
   * `pyyaml`
   * `scipy.optimize`

   Optional packages:
   * `matplotlib.pyplot`
   * `doxygen`

## Usage
```bash
python pd_dic.py -i input.yaml -t pd   
```   
Where `-i` has to be the confriguration in yaml format and `-t` is the type, which can be `pd` for peridynamic simulations.

## Input description

### Material 

The material parameters are described here
```yaml
Material:
    Type: Elastic
    E_Modulus: 4000.0
```
The available `Type` are until now `Elastic` and `Viscoelastic`. 
### Geometry
The discretization and the nodes are described with
```yaml
Discretization:
    Dim: 1
    Final_Time: 2.0
    Time_Steps: 8
    Horizon_Factor_m_value: 1.0
    Influence_Function: 1.0
    File: 
        Name: geometry_dx0_50.csv
```
where `Dim` is the dimension of the node cloud, `Final_Time` the end time of the simulation, `Time_Steps` the amount of time steps, `Horizon_Factor_m_value` the m value of the horizon, `Influence_Function` the factor to scale the influence of the force withrespect to the distance of the horizon, and `Name` the file providing the node information in the CSV format with spaces as delimiter. An example for this file is provided here:
```yaml
#id x y z volume
0 0.0 0.0 0.0 1.0 
1 1.0 1.0 1.0 1.0
```
The shape for the load is given here
```yaml
 Shape:
        Type: Ramp
        Values:
            - 1.5
            - 2.0
            - 2.0
```yaml
where `Type` describes the shape and `Values` specifiy the geomerty of the shape.

### Boundary Conditions
Boundary conditions can be described with
```yaml
Boundary:
    Condition:
        Type: 
            - Force or Displacement
        Value: 
            - Float
        File: 
            - file.csv
```
where the `Type` either can be `Force` or `Displacement`, `Value` describes the value in Newton or Millimeter whis is applied at the nodes
described in `File`. The file has to be provided in the CSV format with spaces as delimiter with the id ofthe nodes where the condition 
should be applied. Here, is an example for a `file.csv`
```yaml
#id
0
1
```
### Output

#### CSV
For writing simulation attributes as CSV format the `Output` tag can be used.
```yaml
Output:
    CSV:
        Type:
            - Position
        File:
            - nodes_positions.csv
```
Where `Type` specifies the attribute and `File` the file name of the output file.

## Examples

An example for an elastic material and an viscoelastic material is provided in the example folder

# License
The code is licensed under the MIT License developed by [Ilyass Tabiai](http://iltabiai.github.io/), Rolland Delorme, and [Patrick Diehl](http://diehlpk.github.io/cv).

Based on a work at <a xmlns:dct="http://purl.org/dc/terms/" href="http://dx.doi.org/10.1016/S0022-5096(99)00029-0" rel="dct:source">http://dx.doi.org/10.1016/S0022-5096(99)00029-0</a>.
