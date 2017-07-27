import numpy as np
from ..util import linalgebra

def w(problem,X,type):
    
    if type == "ONE":
        return 1.
    if type == "EXP":
        len =  linalgebra.norm(X)
        return np.exp(- (len*len) / problem.neighbors.horizon)
    
    return 1.