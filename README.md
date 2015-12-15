# Getting started

###Dependencies

   The following python packages are **required**:
   * `numpy`
   * `xmltodict`
   * `scipy.optimize`

   Optional packages:
   * `matplotlib.pyplot`

###Basic script
Solve a default problem provided in the `deck.xml` by either executing the following python script in the folder where you cloned the repo, or opening a python console in the same folder and executing the following lines one by one.

```python
import numpy as np
import random

#Load the PD_deck class and create a PD_deck object
from deck import PD_deck
data = PD_deck()

#Load the PD_problem class and create a PD_problem object
from problem import PD_problem
problem = PD_problem( data )

#Create an initial guess vector here based on a linear distribtuion of the 
#nodes, disturbed by a small random coefficient
y = np.zeros( ( int(data.Num_Nodes) ) )
#data.Num_Nodes is the total number of nodes provided by the user
for x_i in range(0, int(data.Num_Nodes)):
    # problem.x contains a vector of linearly distributed nodes on the bar
    y[x_i] = problem.x[x_i]+0.1*random.random()*data.Delta_x
#y is our initial guess

#Load the elastic_material class and compute first step PD forces
from elastic import elastic_material
forces = elastic_material( data, problem, y )

#Solve the problem
problem.quasi_static_solver( y, data, forces )

#Write the results to a CSV file
problem.write_data_to_csv(data, problem)
#The problem resolution (time step by time step) is now written in a 
#csv file called data_csv in the current folder

```


## PD_deck class

Open a python terminal or create a python script in the folder where the `deck.py` file is:

   `from deck import PD_deck` will import the *PD_deck* class.

   `data = PD_deck()` will create an object called `data` having the variables and methods of the *PD_deck* class. To create that object, the class will import data from the `deck.xml` file. It is necessary to properly fill it. The following sections explain how to do so.

It is possible to check the data currently loaded in the class using, for example, `data.Final_time` will show you the total duration of the simulation, or `data.Num_Nodes` to see how many PD nodes there are.

#### PD_deck variables

   List of available parameters of the *PD_deck* class:

   * `PD_deck.Horizon_Factor:` Horizon = Horizon\_Factor X `PD_deck.Delta_x` X `Safety parameter`

      * `Safety parameter`= 1.01

   * `PD_deck.N_Delta_t:` Total number of timesteps.

   * `PD_deck.Final_Time:` Total duration of the simulation.

   * `PD_deck.N_Delta_x:` Total number of spatial steps.

   * `PD_deck.Length:` Length of the 1D bar.

   * `PD_deck.Num_Nodes:` Total number od PD nodes.

   * `PD_deck.Delta_x:` Spatial step.

   * `PD_deck.Length_Tot:` Length of the bar specified by the user + additional nodes to apply the load.

   * `PD_deck.x:` Vector x containing the positions of the PD nodes centered around 0 and separated by a constant distance `PD_deck.Delta_x`.

   * `PD_deck.Influence_Function:` Value fo the influence function provided in `deck.xml`

   * `PD_deck.Loading_Flag:` Loading scheme selected by the user.

   * `PD_deck.Material_Flag:` Material behaviour selected by the user.

## PD_problem class

   `from problem import PD_problem` will import the *PD_problem* class.
   
   `problem = PD_problem( data )` will create a `problem` object which is a *PD_prolem* class.
   
   The problem is now loaded and it is possible to check some of its parameters 
   ( `PD_problem.b`, `PD_problem.Horizon`, `PD_problem.x` ). In order to solve the problem, it is necessary to select a material behavior class and provide an initial guess vector.
   
#### PD_problem variables

   * `PD_problem.b` External load vector **b**. Its length is equal to the total number of PD nodes `PD_deck.Num_Nodes`

   * `PD_problem.Horizon` This parameter is Horizon = Horizon\_Factor X `PD_deck.Delta_x` X `Safety parameter`
    * `Safety parameter`= 1.01
    
   * `PD_problem.x` A vector **x** of linearly distributed points on the bar. Its length is equal to `PD_deck.Num_Nodes`

In order to **solve the problem** it is necessary to select a material parameter and load it in an object called `forces` for example. This is covered in the section [elastic_material](https://github.com/joydisee/peridynamics_viscoelasticity_1D#elastic_material-class) class. Solving the problem is covered in [PD_problem methods](https://github.com/joydisee/peridynamics_viscoelasticity_1D#pd_problem-methods) section.

**After solving the problem**, the following variables are also available:

   * `PD_problem.y` A matrix of size `PD_deck.Num_Nodes`x`PD_deck.Num_TimeStep` which records the position of each PD node at each time step.

   * `PD_problem.forces` A matrix of size `PD_deck.Num_Nodes`x`PD_deck.Num_TimeStep` which records the PD forces at each PD node, at each time step.

#### PD_problem methods

   * `PD_problem.quasi_static_solver( y, PD_deck, elastic_material )` Takes an initial guess vector **y** of length `PD_deck.Num_Nodes`, a *PD_deck* class object and an *elastic_material* class ibject. **It solves the problem** and provides the `PD_problem.y` and `PD_problem.forces` variables.
   
   * `PD_problem.write_data_to_csv( PD_deck, PD_problem )` Takes a *PD_deck* class object and a *PD_problem* class object and writes the result of the solved problem in a csv file `data_csv`. This method is **not available before solving the problem**.

## elastic_material class



## Example of XML deck

```XML
<?xml version="1.0"?>
<data>
	<Discretization>
			<N_Delta_x>4</N_Delta_x>
			<N_Delta_t>2</N_Delta_t>
			<Final_Time>1.0</Final_Time>
			<Horizon_Factor>1</Horizon_Factor>
			<Influence_Function>1</Influence_Function>
	</Discretization>
	<Boundary_Conditions>
	        <Type>RAMP</Type>
			<Force>54.0</Force>
			<Ramp_Time>1.0</Ramp_Time>
	</Boundary_Conditions>
	<Geometry>
			<Length>4.0</Length>
			<Surface>1</Surface>
	</Geometry>
	<Material>
	        <Type>ELASTIC</Type>
			<E_Modulus>4000.0</E_Modulus>
	</Material>
</data>
```

`<Discretization>` must contain:

* `<N_Delta_x>` The number of spatial steps.
* `<N_Delta_t>` The number of time steps.
* `<Final_Time>` Duration of the simulation.
* `<Influence_Function>` The value of the influence function as defined in Silling, Lehoucq 2010.

`<Boundary_Conditions>` must contain:

* `<Type>` Currently, only the `RAMP` loading type is available.
    * `<Force>` Force applied on the bar on its extremeties.
    * `<Ramp_Time>` Time during which the force grows linearly before reaching the force value.

`<Geometry>` must contain:

* `<Length>` Lenth of the 1D bar.
* `<Surface>` Surface occupied by the 1D bar (perpendicular to the bar's direction).

`<Material>` must contain:

* `<Type>` Currently, only the `ELASTIC` material is avilable.
    * `<E_Modulus>` Value of the 1D elastic modulus.
*
