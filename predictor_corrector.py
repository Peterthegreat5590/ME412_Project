
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
        uNW = u[2:,:-2]

        vP = v[1:-1,1:-1]
        vE = v[1:-1,2:]
        vW = v[1:-1,:-2]
        vN = v[2:,1:-1]
        vS = v[:-2,1:-1]
        vSE = v[:-2,2:]

        ue = (uE+uP)/2
        uw = (uW+uP)/2

        vn_u = (vP+vE)/2
        vs_u = (vS+vSE)/2

        vn = (vN+vP)/2
        vs = (vS+vP)/2

        ue_v = (uP+uN)/2
        uw_v = (uW+uNW)/2


        Re_u = np.abs(uP)*dy/nu
        Re_v = np.abs(vP)*dx/nu


        a_uE = nu + np.maximum(-ue*dy,0)
        a_uW = nu + np.maximum(uw*dy,0)
        a_uN = nu + np.maximum(-vn_u*dx,0)
        a_uS = nu + np.maximum(vs_u*dx,0)
        a_uP = a_uE + a_uW + a_uN + a_uS + (ue*dy - uw*dy) + (vn_u*dx - vs_u*dx)

        a_vE = nu + np.maximum(-ue_v*dy,0)
        a_vW = nu + np.maximum(uw_v*dy,0)
        a_vN = nu + np.maximum(-vn*dx,0)
        a_vS = nu + np.maximum(vs*dx,0)
        a_vP = a_vE + a_vW + a_vN + a_vS + (ue_v*dy - uw_v*dy) + (vn*dx - vs*dx)


        a_uE[Re_u<=2] = nu - ue[Re_u<=2]*dy/2
        a_uW[Re_u<=2] = nu + uw[Re_u<=2]*dy/2
        a_uN[Re_u<=2] = nu - vn_u[Re_u<=2]*dx/2
        a_uS[Re_u<=2] = nu + vs_u[Re_u<=2]*dx/2
        a_uP[Re_u<=2] = a_uE[Re_u<=2] + a_uW[Re_u<=2] + a_uN[Re_u<=2] + a_uS[Re_u<=2] + (ue[Re_u<=2]*dy - uw[Re_u<=2]*dy)/2 + (vn_u[Re_u<=2]*dx - vs_u[Re_u<=2]*dx)/2

        a_vE[Re_v<=2] = nu - ue_v[Re_v<=2]*dy/2
        a_vW[Re_v<=2] = nu + uw_v[Re_v<=2]*dy/2
        a_vN[Re_v<=2] = nu - vn[Re_v<=2]*dx/2
        a_vS[Re_v<=2] = nu + vs[Re_v<=2]*dx/2
        a_vP[Re_v<=2] = a_vE[Re_v<=2] + a_vW[Re_v<=2] + a_vN[Re_v<=2] + a_vS[Re_v<=2] + (ue_v[Re_v<=2]*dy - uw_v[Re_v<=2]*dy)/2 + (vn[Re_v<=2]*dx - vs[Re_v<=2]*dx)/2



        u_star[1:-1,1:-1] = (a_uE*uE + a_uW*uW + a_uN*uN + a_uS*uS)/a_uP
        v_star[1:-1,1:-1] = (a_vE*vE + a_vW*vW + a_vN*vN + a_vS*vS)/a_vP

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

        

        if t%(1000*dt)<dt/2:
            np.savetxt(f"u_Time={t}.csv", u, '%0.5f', ',','\n')
            np.savetxt(f"v_Time={t}.csv", v, '%0.5f', ',','\n')
            np.savetxt(f"p_Time={t}.csv", p, '%0.5f', ',','\n')
        
        if t%(dt*100)<dt/2:
            print(f"Time={t}")

    np.savetxt(f"u_Time={t}.csv", u, '%0.5f', ',','\n')
    np.savetxt(f"v_Time={t}.csv", v, '%0.5f', ',','\n')
    np.savetxt(f"p_Time={t}.csv", p, '%0.5f', ',','\n')


