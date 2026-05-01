
import numpy as np
from sor_solve import sor_solve
from enforce_boundary import enforce_boundary
from SOR_solve_copy import SORitersolve


def solve(u, v, p, u_mask, v_mask, p_mask, dt, dx, dy, nu, X_p, Y_p, X_u, Y_u, X_v, Y_v, L, D, Vin, time):

    u_star = np.zeros_like(u)
    v_star = np.zeros_like(v)

  
    for t in time:

        uP = u[1:-1,1:-1]
        uE = u[1:-1,2:]
        uW = u[1:-1,:-2]
        uN = u[2:,1:-1]
        uS = u[:-2,1:-1]

        ue = (uE+uP)/2
        uw = (uW+uP)/2
        un = (uN+uP)/2
        us = (uS+uP)/2

        vn = (v[1:-1,1:-1]+v[1:-1,2:])/2
        vs = (v[:-2,1:-1]+v[:-2,2:])/2

        u_conv = np.maximum(-ue,0)*dy*ue + np.maximum(uw,0)*dy*uw + np.maximum(-vn,0)*dx*un + np.maximum(vs,0)*dx*us #Negative already applied by signs, Upwinding Scheme Used

        u_diff = nu*(uE+uW+uN+uS-4*uP)

        u_star[1:-1,1:-1] = uP + dt/(dx*dy)*(u_conv+u_diff)

        vP = v[1:-1,1:-1]
        vE = v[1:-1,2:]
        vW = v[1:-1,:-2]
        vN = v[2:,1:-1]
        vS = v[:-2,1:-1]

        ve = (vE+vP)/2
        vw = (vW+vP)/2
        vn = (vN+vP)/2
        vs = (vS+vP)/2

        ue = (u[1:-1,1:-1]+u[2:,1:-1])/2
        uw = (u[1:-1,:-2]+u[2:,:-2])/2

        v_conv = np.maximum(-ue,0)*dy*ve + np.maximum(uw,0)*dy*vw + np.maximum(-vn,0)*dx*vn + np.maximum(vs,0)*dx*vs #Negative already applied by signs, Upwinding Scheme Used

        v_diff = nu*(vE+vW+vN+vS-4*vP)

        v_star[1:-1,1:-1] = vP + dt/(dx*dy)*(v_conv+v_diff)

        u_star, v_star, _ = enforce_boundary(u_star, v_star, p, u_mask, v_mask, p_mask, X_p, Y_p, X_u, Y_u, X_v, Y_v, L, D, Vin)

        p = sor_solve(u_star, v_star, p, dx, dy, dt, X_p, Y_p, L, D)

        u_new = u_star.copy()
        u_new[1:-1,1:-1] = u_star[1:-1,1:-1] - dt*(p[1:-1,2:]-p[1:-1,1:-1])/dx

        v_new = v_star.copy()
        v_new[1:-1,1:-1] = v_star[1:-1,1:-1] - dt*(p[2:,1:-1]-p[1:-1,1:-1])/dy

        u[u_mask] = u_new[u_mask]
        v[v_mask] = v_new[v_mask]

        u, v, p = enforce_boundary(u, v, p, u_mask, v_mask, p_mask, X_p, Y_p, X_u, Y_u, X_v, Y_v, L, D, Vin)

        if np.any(np.isnan(u)):
            print(t)

        

        if np.isclose(t%(100*dt),0):
            np.savetxt(f"u_Time={t}.csv", u, '%0.5f', ',','\n')
            np.savetxt(f"v_Time={t}.csv", v, '%0.5f', ',','\n')
            np.savetxt(f"p_Time={t}.csv", p, '%0.5f', ',','\n')



