from numba import *
import numpy as np

@jit
def SORitersolve(p, ux, uy, rho, dt, dx, dy, omega, tol, maxiter) :

    f = np.zeros_like(p)
    f[1:-1,1:-1] = dx*dy*rho/dt*((ux[1:-1,2:]-ux[1:-1,:-2])/(2*dx)+(uy[2:,1:-1]-uy[:-2,1:-1])/(2*dy))
    Nx = p.shape[1]-1
    Ny = p.shape[0]-1
    # Track convergence
    dvals = np.zeros(maxiter)

    converged = False
    iter = 0
    while(not converged):

        # The 'old' guess (iteration k) is saved as the result from the previous guess
        p_old = np.copy(p)

        # Iteration counter
        iter += 1

        # Loop over all interior nodes
        for i in range(1, Ny):
            for j in range(1, Nx):

                # # Interior nodes affected by the Neumann condition: i=N-1
                p[i,j] = p[i,j] + omega*((1/4)*(p[i+1,j] + p[i-1,j] + p[i,j+1] + p[i,j-1] - f[i,j]) - p[i,j])

        p = p - np.mean(p)
        # Check if converged (apply norm to interior nodes only)
        d = np.linalg.norm(p[1:-1, 1:-1] - p_old[1:-1, 1:-1])
        dvals[iter-1] = d

        # Accept or reject current solution
        if(d < tol):
            converged = True
            # Impose the Neumann condition on the edges
            p[Ny,:] = p[Ny-1,:]
            p[0,:] = p[1,:]
            p[:,Nx] = p[:,Nx-1]
            p[:,0] = p[:,1]


        elif(iter==maxiter):
            print('no convergence')
            break
    
    dvals = dvals[dvals!=0]

    return p, dvals