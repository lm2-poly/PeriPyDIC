from math import sqrt 

def norm(a):
    
    if len(a.shape) == 1:
        if a.shape[0] ==  1:
            return abs(a)
    
        if a.shape[0] ==  2:
            return sqrt(a[0]*a[0]+a[1]*a[1])
        
    if len(a.shape) == 2:
        if a.shape[1] ==  2:
            tmp = 0. 
            for i in range(0,a.shape[0]):
                tmp += a[i][0]*a[i][0]+a[i][1]*a[i][1]
            return sqrt(tmp) 
        if a.shape[1] == 1:
            return sum(abs(a))