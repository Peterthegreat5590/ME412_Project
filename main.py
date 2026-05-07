import numpy as np
import sys
import matplotlib.pyplot as plt
from enforce_boundary import enforce_boundary
from predictor_corrector import solve

inputs = sys.argv

if len(inputs)>4:
    raise ValueError("Invalid input")
elif len(inputs)<4:
    dx = 1/16
    dy = 1/16
    dt = 1/10000
    t_max = 5
elif "/" in sys.argv:
    dx = float(eval(sys.argv[1]))
    dy = dx
    dt = float(sys.argv[2])
    t_max = float(sys.argv[3])
else:
    raise ValueError("Invalid input")

L = 5 #Length of channel
D = 1 #Witdh of channel

H = L+D
W = L
Vin = 50
Re=50
nu = Vin*D/Re


x_p = np.arange(-dx, W+dx*3/2, dx) + dx/2
y_p = np.arange(-dy, H+dy*3/2, dy) + dy/2
x_u = np.arange(-dx, W+dx*3/2, dx) + dx
y_u = np.arange(-dy, H+dy*3/2, dy) + dy/2
x_v = np.arange(-dx, W+dx*3/2, dx) + dx/2
y_v = np.arange(-dy, H+dy*3/2, dy) + dy

Nx = x_p.shape[0]-2
Ny = y_p.shape[0]-2

u = np.zeros([Ny+2, Nx+2])
v = np.zeros([Ny+2, Nx+2])
p = np.zeros([Ny+2, Nx+2])

X_p, Y_p = np.meshgrid(x_p,y_p)
X_u, Y_u = np.meshgrid(x_u,y_u)
X_v, Y_v = np.meshgrid(x_v,y_v)

domain_p = np.ones([Ny+2, Nx+2],bool)
domain_p[(X_p>D)&(X_p<(W-D))&(Y_p>D)] = False


domain_u = np.ones([Ny+2, Nx+2],bool)
domain_u[((X_u>D)&(X_u<(W-D))&(Y_u>D))|(X_u>W+dx)] = False

domain_v = np.ones([Ny+2, Nx+2],bool)
domain_v[((X_v>D)&(X_v<(W-D))&(Y_v>D))|(Y_v>H+dy)] = False


time = np.arange(0,t_max,dt)

u, v, p = enforce_boundary(u, v, p, domain_u, domain_v, domain_p, X_p, Y_p, X_u, Y_u, X_v, Y_v, L, D, Vin)

np.set_printoptions(threshold=100000)

np.savetxt("X_p.csv",X_p,'%0.5f',',','\n')
np.savetxt("Y_p.csv",Y_p,'%0.5f',',','\n')
np.savetxt("X_u.csv",X_u,'%0.5f',',','\n')
np.savetxt("Y_u.csv",Y_u,'%0.5f',',','\n')
np.savetxt("X_v.csv",X_v,'%0.5f',',','\n')
np.savetxt("Y_v.csv",Y_v,'%0.5f',',','\n')

solve(u, v, p, domain_u, domain_v, domain_p, dt, dx, dy, nu, X_p, Y_p, X_u, Y_u, X_v, Y_v, L, D, Vin, time)

