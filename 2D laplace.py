import numpy as np
import matplotlib.pyplot as plt


#Laplace 2-D Function

def laplace2d(p, y, dx, dy, l1norm_target):
    l1norm = 1
    pn = np.empty_like(p)

    while l1norm > l1norm_target:
        pn = p.copy() #Create a new list and saves old list if while condition is not satisfied
        p[1:-1,1:-1] = ((dy ** 2 * (pn[1:-1, 2:] + pn[1:-1, 0:-2]) + dx ** 2 * (pn[2:, 1:-1] + pn[0:-2, 1:-1])) / (2 *
                                                        (dx ** 2 + dy ** 2))) #Discretized equation to solve for p

        #Boundary Conditions
        p[:, 0] = 0 #p = 0 @ x = 0
        p[:,-1] = y #p = y @ x = 2
        p[0, :] = p[1, :] # dp/dy = 0 @ y = 0
        p[-1,:] = p[-2, :] # dp/dy = 0 @ y = 1
        l1norm = (np.sum(np.abs(p[:]) - np.abs(pn[:]))) / (np.sum(np.abs(pn[:]))) #Compares copy of old p 'pn' to new p
    return p



