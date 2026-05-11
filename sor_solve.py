from numba import *
import numpy as np
from enforce_boundary import pressure_boundary

@jit
def sor_solve(u, v, p, dx, dy, dt, X_p, Y_p, L, D, p_mask):
    omega = 1.8
    b = dx/dt*(u[1:-1,1:-1]-u[1:-1,:-2]+v[1:-1,1:-1]-v[:-2,1:-1])
    converged = False
    tolerance = 1e-5
    max_iters = 10000
    
    for iter in range(max_iters):
        p_old = p.copy()
        for i in range(1,u.shape[0]-1):
            for j in range(1,u.shape[1]-1):
                if not p_mask[i,j]:
                    continue
                p_new = ((p[i,j+1] + p[i,j-1] + p[i+1,j] + p[i-1,j]) - b[i-1,j-1])/4
                p[i,j] = p[i,j] + (p_new - p[i,j])*omega
        diff = np.max(np.abs(p_old[1:-1,1:-1]-p[1:-1,1:-1]))
        #print(p)
        converged = diff<tolerance
        if converged:
            p = pressure_boundary(p, X_p, Y_p, L, D)
            break
            
    return p
        

        
    