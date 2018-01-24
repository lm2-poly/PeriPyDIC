import numpy as np
from ..util import linalgebra

def w(problem,X,type):
    
    if type == "ONE":
        return 1.
    if type == "EXP":
        len =  linalgebra.norm(X)
        return np.exp(- (len*len) / problem.neighbors.horizon / problem.neighbors.horizon)
    if type == "NORM":
        return 1. / linalgebra.norm(X) 
    
    return 1.

def damage(deck,problem,i,j):
    
    if deck.damage_type == "ADAPTIVE":
        return (1-problem.neighbors.damage[i][j])
    
    if deck.damage_type == "FULL":
        
        if (problem.neighbors.damage[i][j] >= 1.):
            return 0.
        else:
            return 1.
        
    return 1.