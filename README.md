# Basic Usage

```python
from peripydic import *

deck = PD_deck("./input_elas_1D.yaml")
problem = PD_problem(deck)

```

# The problem

   Consider a 1D bar. The length of the bar and the number of sections in which the bar will be divided (the Peridynamic (PD) nodes are at the edges of each section) are provided by the user.

   A ramp loading is applied at the edges of the bar. The user provides the number of seconds after which the ramp reaches a maximal force value and remains there.
   Additional blocks are added at the edges of the block in order to apply the load. The number of blocks added at each edge of the bar is provided as the Horizon Factor.
   The user also provides the total duration of the simulation and the total number of time steps desired.

   The problem is then solved using an initial guess vector and the *Newton-Krylov* algorithm provided in the `scipy.optimize` package.  

# Getting started

### Dependencies

   The following python packages are **required**:
   * `numpy`
   * `pyyaml`
   * `scipy.optimize`
   * `sharedmem`

   The following tools are **optional**
   * `doxygen`
   * `dot`

### Installation

    
```bash
virtualenv pddic
source pddic/binactivate
pip install -r requirements.txt
python setup.py install
```

## Usage
```bash
python pd_dic.py -i input.yaml -t pd   
```   
Where `-i` has to be the configuration in `yaml format` and `-t` is the type, which can be `pd` for peridynamic simulations and `dic` for processing results from digital image correlation.

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
    Saftety_Factor: 1.001
    File:
        Name: geometry_dx0_50.csv
        Path: ./
        Type: mudic
```
where `Dim` is the dimension of the node cloud, `Final_Time` the end time of the simulation, `Time_Steps` the amount of time steps, `Horizon_Factor_m_value` the m value of the horizon, `Influence_Function` the factor to scale the influence of the force with respect to the distance of the horizon and `Saftety_Factor` influences the computation of the horizon, and `Name` the file providing the node information in the CSV format with spaces as delimiter. The `path` is the path to the file in the file system. The `Type` describes if the CSV file was exported by `mudic` or `vic3d`. An example for this file is provided here:

```
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
```
where `Type` describes the shape and `Values` specify the geometry of the shape.

### Boundary Conditions

Boundary conditions can be described with
```yaml
Boundary:
    Condition:
        Type:
            - Force or Displacement
        Value:
            - Float
        Direction:
            - Int
        File:
            - file.csv
```
where the `Type` either can be `Force` or `Displacement`, `Value` describes the value in Newton or Millimeter whis is applied at the nodes
described in `File`, and `Direction` describes the direction (X=1,Y=2,Z=3) where the condition is applied. The file has to be provided in the CSV format with spaces as delimiter with the id of the nodes where the condition
should be applied. Here, is an example for a `file.csv`
```yaml
#id
0
1
```
### Output

For writing simulation attributes the `Output` tag can be used.

#### CSV

For writing the simulation attributes to the CSV format the tag `CSV` is used.

```yaml
Output:
    CSV:
        Type:
            - Position
        File:
            - nodes_positions.csv
```
Where `Type` specifies the attribute and `File` the file name of the output file.

#### VTK

For writing the simulation attributes to the VTK unstructured grid format the tag `VTK` is used

```yaml
VTK:
    Path: ./
    Type:
    - Displacement
    - Neighbors
    - Force
    - Conditions
    - Volume_Force
    - Strain
Slice: 1
```
Where `Path` is the path for the output, `Type` specify the simulation attributes, which are considered for the output, and `Slice` defines that every n-th time step is written.

### Solver

Here, `Max_Iteration`, `Tolerance` of the solver can be specified. With `Jacobian_Perturbation` the perturbation for assembly the Jacobian matrix is defined.

```yaml
Solver:
    Max_Iteration: 100
    Tolerance: 1.0e-6
    Jacobian_Perturbation: 1.0e-6
```

### Parallel computing
For using multiple threads with `multiprocessing` specify the number of threads with `Threads`.

```yaml
Parallel:
    Threads: 2
```


## Examples

An example for an elastic material and an viscoelastic material is provided in the example folder

# Publications

* Delorme, R., Tabiai, I., Laberge Lebel, L., & Lévesque, M. (2017). **Generalization of the ordinary state-based peridynamic model for isotropic linear viscoelasticity. Mechanics of Time-Dependent Materials.**, _Mechanics of Time-Dependent Materials, 1-27_, 10.1007/s11043-017-9342-3,
* Rolland Delorme, Patrick Diehl, Ilyass Tabiai, Louis Laberge Lebel, and Martin
Lévesque. **Extracting constitutive mechanical parameters in linear elasticity using the virtual fields method within the ordinary state-based peridynamic framework**. Journal of Peridynamics and Nonlocal Modeling, Jan 2020. [Link](https://link.springer.com/article/10.1007%2Fs42102-019-00025-7), [Preprint](https://engrxiv.org/uv8m7/)

# License

The code is licensed under the GNU General Public License v3.0 developed by [Patrick Diehl](http://diehlpk.github.io/), [Rolland Delorme](https://orcid.org/0000-0001-7637-3936) and [Ilyass Tabiai](http://iltabiai.github.io/) . Please cite our code with following [![DOI](https://zenodo.org/badge/46075533.svg)](https://zenodo.org/badge/latestdoi/46075533)

Based on works at <a xmlns:dct="http://purl.org/dc/terms/" href="http://dx.doi.org/10.1016/S0022-5096(99)00029-0" rel="dct:source">http://dx.doi.org/10.1016/S0022-5096(99)00029-0</a> and <a xmlns:dct="http://purl.org/dc/terms/" href="https://doi.org/10.1007/s11043-017-9342-3" rel="dct:source">https://doi.org/10.1007/s11043-017-9342-3</a>.
