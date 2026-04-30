def enforce_boundary(u, v, p, u_mask, v_mask, p_mask, X_p, Y_p, X_u, Y_u, X_v, Y_v, L, D, Vin):


    u[not u_mask] = 0
    v[not v_mask] = 0
    p[not p_mask] = 0
    dx = X_p[0,1] - X_p[0,0]
    dy = Y_p[1,0] - Y_p[0,0]


    # Side Boundaries
    u[(X_u==0)|(X_u==L)] = 0
    v[(X_v==-dx/2)] = -v[(X_v==dx/2)]
    v[(X_v==L+dx/2)] = -v[(X_v==L-dx/2)]
    
    # Bottom Boundary
    u[(Y_u==-dy/2)] = -u[(Y_u==dy/2)]
    v[(Y_v==0)] = 0

    # Central Side Boundaries
    u[(X_u==D)|(X_u==L-D)] = 0
    v[(X_v==L-D-dx/2)] = -v[(X_v==L-D+dx/2)]
    v[(X_v==D+dx/2)] = -v[(X_v==D-dx/2)]

    # Central Bottom Boundary
    u[(Y_u==D+dy/2)] = -u[(Y_u==D-dy/2)]
    v[(Y_v==D)] = 0

    # Inlet/Outlet Boundary
    u[(Y_u==D+L+dy/2)&(X_v>0)&(X_v<D)] = -u[(Y_u==D+L-dy/2)&(X_v>0)&(X_v<D)]
    v[(Y_v==D+L)&(X_v>0)&(X_v<D)] = -Vin

    v[(Y_v==D+L)&(X_v>L-D)&(X_v<L)] = v[(Y_v==D+L-dy)&(X_v>L-D)&(X_v<L)]
    u[(Y_u==D+L+dy/2)&(X_v>L-D)&(X_v<L)] = u[(Y_u==D+L-dy/2)&(X_v>L-D)&(X_v<L)]

    p[(Y_p==D+L+dy/2)&(X_p>L-D)&(X_p<L)] = -p[(Y_p==D+L-dy/2)&(X_p>L-D)&(X_p<L)]


    return u, v, p