# -*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrick.diehl@polymtl.ca
import util.neighbor
from materials.elastic import Elastic_material
import numpy as np

## Copmutes extension and force states out of displacement/nodes data obtained
# from the VIC3D CSV (DIC data)
class DIC_problem():

    ## Constructor
    # Find neighbors for DIC data, computes the weighted function, then computes
    # actual positions, extension states and force states
    # @param deck Deck object containig data from the .yaml file
    def __init__(self, deck):
        ## NeighborSearch
        self.neighbors = util.neighbor.NeighborSearch(deck)


        ## Compute the weighted volume for each node in a vector.
        self.weighted_function(deck)

        ## Actual position from DIC result
        self.y = np.zeros((deck.num_nodes, deck.dim,1),dtype=np.float64)

        ## Internal forces
        self.force_int = np.zeros((deck.num_nodes, deck.dim,1),dtype=np.float64)
        ## Extension state
        self.ext = np.zeros( ( deck.num_nodes, deck.num_nodes,1),dtype=np.float64 )


        if deck.material_type == "Elastic":

            mat_class = Elastic_material( deck, self, deck.geometry.act )
            self.update_force_data(mat_class)
            self.update_ext_state_data(mat_class)

        self.update_pos(deck.geometry.act)

    ## Computes the weights for each PD node
    # @param deck Deck object containig data from the .yaml file
    def weighted_function(self, deck):
        ## Weighted volumes vector
        self.weighted_volume = np.zeros((deck.num_nodes),dtype=np.float64)
        for i in range(0, deck.num_nodes):
            index_x_family = self.neighbors.get_index_x_family(i)
            for p in index_x_family:
                X = deck.geometry.nodes[p,:] - deck.geometry.nodes[i,:]
                self.weighted_volume[i] += deck.influence_function * (np.linalg.norm(X))**2 * deck.geometry.volumes[p]

    ## Records the force vector at each time step
    # @param mat_class Material class object for the elastic/viscoelastic material models
    def update_force_data(self, mat_class):
        ## Internal forces state
        self.force_int[:,:, 0] = mat_class.f_int


    ## Records the ext_state vector at each time step
    # @param mat_class Material class object for the elastic/viscoelastic material models
    def update_ext_state_data(self, mat_class):
        ## Extension state
        self.ext[:, :, 0] = mat_class.e

    ## Records the actual position vector at each time step
    # @param act Actual position obtained from DIC data
    def update_pos(self,act):
        ## Actual position state
        self.y[:,:, 0] = act
