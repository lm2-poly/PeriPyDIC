import numpy as np

def w(deck,X,type):
    
    if type == "ONE":
        return 1.
    if type == "EXP":
        len =  X.length()
        return np.exp(- (len*len) / deck.geometry.neighbors.horizon)
    
    return 1.