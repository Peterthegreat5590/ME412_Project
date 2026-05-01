from numba import *
import numpy as np

def enforce_boundary(u, v, p, u_mask, v_mask, p_mask, X_p, Y_p, X_u, Y_u, X_v, Y_v, L, D, Vin):


    u[np.logical_not(u_mask)] = 0
    v[np.logical_not(v_mask)] = 0

    dx = X_p[0,1] - X_p[0,0]
    dy = Y_p[1,0] - Y_p[0,0]


    # Side Boundaries
    u[(X_u==0)|(X_u==L)] = 0
    v[(X_v==-dx/2)] = -v[(X_v==dx/2)]
    v[(X_v==L+dx/2)] = -v[(X_v==L-dx/2)]
    p[(X_p==L+dx/2)] = p[(X_p==L-dx/2)]
    p[(X_p==-dx/2)] = p[(X_p==dx/2)]
    
    # Bottom Boundary
    u[(Y_u==-dy/2)] = -u[(Y_u==dy/2)]
    v[(Y_v==0)] = 0
    p[(Y_p==-dy/2)] = p[(Y_p==dy/2)]

    # Central Side Boundaries
    u[((X_u==D)|(X_u==L-D))&(Y_u>D)] = 0
    v[(X_v==L-D-dx/2)&(Y_v>D)] = -v[(X_v==L-D+dx/2)&(Y_v>D)]
    v[(X_v==D+dx/2)&(Y_v>D)] = -v[(X_v==D-dx/2)&(Y_v>D)]
    p[(X_p==L-D-dx/2)&(Y_p>D)] = p[(X_p==L-D+dx/2)&(Y_p>D)]
    p[(X_p==D+dx/2)&(Y_p>D)] = p[(X_p==D-dx/2)&(Y_p>D)]

    # Central Bottom Boundary
    u[(Y_u==D+dy/2)&(X_u>D)&(X_u<L-D)] = -u[(Y_u==D-dy/2)&(X_u>D)&(X_u<L-D)]
    v[(Y_v==D)&(X_v>D)&(X_v<L-D)] = 0
    p[(Y_p==D+dy/2)&(X_p>D)&(X_p<L-D)] = p[(Y_p==D-dy/2)&(X_p>D)&(X_p<L-D)]

    # Inlet/Outlet Boundary
    u[(Y_u==D+L+dy/2)&(X_v>0)&(X_v<D)] = -u[(Y_u==D+L-dy/2)&(X_v>0)&(X_v<D)]
    v[(Y_v==D+L)&(X_v>0)&(X_v<D)] = -Vin

    v[(Y_v==D+L)&(X_v>L-D)&(X_v<L)] = v[(Y_v==D+L-dy)&(X_v>L-D)&(X_v<L)]
    u[(Y_u==D+L+dy/2)&(X_v>L-D)&(X_v<L)] = u[(Y_u==D+L-dy/2)&(X_v>L-D)&(X_v<L)]

    p[(Y_p==D+L+dy/2)&(X_p>L-D)&(X_p<L)] = -p[(Y_p==D+L-dy/2)&(X_p>L-D)&(X_p<L)]


    return u, v, p

@jit
def pressure_boundary(p, X_p, Y_p, L, D):
    dx = X_p[0,1] - X_p[0,0]
    dy = Y_p[1,0] - Y_p[0,0]

    # Outlet Boundary
    for i in range(p.shape[0]):
        for j in range(p.shape[1]):
            if (Y_p[i,j]==D+L+dy/2)&(X_p[i,j]>L-D)&(X_p[i,j]<L):
                p[i,j] = -p[i-1,j]

            # Central U Boundary
            if (Y_p[i,j]==D+dy/2)&(X_p[i,j]>D)&(X_p[i,j]<L-D):
                p[i,j] = p[i-1,j]
            if (X_p[i,j]==D+dx/2)&(Y_p[i,j]>D):
                p[i,j] = p[i,j-1]
            if (X_p[i,j]==L-D-dx/2)&(Y_p[i,j]>D):
                p[i,j] = p[i,j+1]
            if (X_p[i,j]==D+dx/2)&(Y_p[i,j]==D+dy/2):
                p[i,j] = (p[i,j-1]+p[i-1,j])/2
            if (X_p[i,j]==L-D-dx/2)&(Y_p[i,j]==D+dy/2):
                p[i,j] = (p[i,j+1]+p[i-1,j])/2

    return p