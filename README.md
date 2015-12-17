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
#nodes, disturbed by a small random coefficient using a method provided in the
#PD_problem class
x_0 = problem.provide_random_initial_guess( data ) 
#x_0 is our initial guess

#Solve the problem, the material behavior model provided in the deck.xml
#will be used
problem.quasi_static_solver( x_0, data )

#Check the position of PD nodes at the 3rd time step
print problem.y[:, 3]

#Check the PD force value at each node at the 5th time step
print problem.forces[:, 5]

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
   ( `PD_problem.b`, `PD_problem.Horizon`, `PD_problem.x` ). In order to solve the problem, it is necessary to provide an initial guess vector.
   
#### PD_problem variables

   * `PD_problem.b` External load vector **b**. Its length is equal to the total number of PD nodes `PD_deck.Num_Nodes`

   * `PD_problem.Horizon` This parameter is Horizon = Horizon\_Factor X `PD_deck.Delta_x` X `Safety parameter`
    * `Safety parameter`= 1.01
    
   * `PD_problem.x` A vector **x** of linearly distributed points on the bar. Its length is equal to `PD_deck.Num_Nodes`

**Solving the problem** is convered in [PD_problem methods](https://github.com/joydisee/peridynamics_viscoelasticity_1D#pd_problem-methods).

**After solving the problem**, the following variables are also available:

   * `PD_problem.y` A matrix of size `PD_deck.Num_Nodes`x`PD_deck.Num_TimeStep` which records the position of each PD node at each time step.

   * `PD_problem.forces` A matrix of size `PD_deck.Num_Nodes`x`PD_deck.Num_TimeStep` which records the PD forces at each PD node, at each time step.

#### PD_problem methods

   * `PD_problem.provide_random_initial_guess( PD_deck ) ` Takes a *PD_deck* class object and **returns** an initial guess vector which is based on the linear distribution of PD nodes on the bar, disturbed by a small random parameter. You can also make your own initial guess vector. It must be of length `PD_deck.Num_Nodes`

   * `PD_problem.quasi_static_solver( x_0, PD_deck )` Takes an initial guess vector **x_0** of length `PD_deck.Num_Nodes` and a *PD_deck* class object. **It solves the problem** and provides the `PD_problem.y` and `PD_problem.forces` variables.
   
   * `PD_problem.write_data_to_csv( PD_deck, PD_problem )` Takes a *PD_deck* class object and a *PD_problem* class object and writes the result of the solved problem in a csv file `data_csv`. This method is **not available before solving the problem**.

## elastic_material class

   `from elastic import elastic_material` will import the *elastic_material* class.
   
   `forces = elastic_material( data, problem, y  )` will create a `forces` object which is an *elastic_material* class. It needs a `data` object of class *PD_deck*, a `problem` object of class *PD_problem* and a position vector **y**. It is not necessary for the user to interact with this class.
   It can be used to compute PD forces for any given positions **y** of the PD nodes. After creating a *PD_deck* object and a *PD_problem* project, you can create a vector **y** of length `PD_deck.Num_Nodes` and then create a `forces = elastic_material( data, problem, y  )`. It is then possible to check the *elastic_material* variables.
   
#### elastic_material class variables

   * `elastic_material.e` A matrix ( `PD_deck.Num_Nodes`, Horizon\_Factor X `PD_deck.Delta_x` ) that provides the deformation of each point relatively to the nodes in its family.
   
   * `elastic_material.T`
   
   * `elastic_material.Ts` A vector of length `PD_deck.Num_Nodes` which contains the PD forces value at each PD node.
   
   

## XML deck description

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

<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/80x15.png" /></a><br /><span xmlns:dct="http://purl.org/dc/terms/" property="dct:title">peridynamics_1D</span> by <a xmlns:cc="http://creativecommons.org/ns#" href="http://iltabiai.github.io/" property="cc:attributionName" rel="cc:attributionURL">Ilyass Tabiai, Rolland Delorme</a> is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.<br />Based on a work at <a xmlns:dct="http://purl.org/dc/terms/" href="http://dx.doi.org/10.1016/S0022-5096(99)00029-0" rel="dct:source">http://dx.doi.org/10.1016/S0022-5096(99)00029-0</a>.
