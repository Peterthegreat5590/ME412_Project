from numba import *
import numpy as np
import sys
import matplotlib.pyplot as plt

inputs = sys.argv

if len(inputs)>2:
    raise("Invalid input")
elif len(inputs)<2:
    dx = 1/64
    dy = 1/64
elif "/" in sys.argv:
    dx = float(eval(sys.argv[1]))
    dy = dx
else:
    raise("Invalid input")

L = 5 #Length of channel
D = 1 #Witdh of channel

H = L+D
W = L

Nx = int(W/dx)
Ny = int(H/dy)

u = np.zeros([Ny+2, Nx+2])
v = np.zeros([Ny+2, Nx+2])
p = np.zeros([Ny+2, Nx+2])
x_p = np.linspace(-dx, W+dx, Nx+2) + dx/2
y_p = np.linspace(-dy, H+dy, Ny+2) + dy/2
x_u = np.linspace(-dx, W+dx, Nx+2) + dx
y_u = np.linspace(-dy, H+dy, Ny+2) + dy/2
x_v = np.linspace(-dx, W+dx, Nx+2) + dx/2
y_v = np.linspace(-dy, H+dy, Ny+2) + dy

X_p, Y_p = np.meshgrid(x_p,y_p)
X_u, Y_u = np.meshgrid(x_u,y_u)
X_v, Y_v = np.meshgrid(x_v,y_v)

domain_p = np.ones([Ny+2, Nx+2])
domain_p[(X_p>D)&(X_p<(W-D))&(Y_p>D)] = 0


domain_u = np.ones([Ny+2, Nx+2])
domain_u[(X_u>D)&(X_u<(W-D))&(Y_u>D)] = 0

domain_v = np.ones([Ny+2, Nx+2])
domain_v[(X_v>D)&(X_v<(W-D))&(Y_v>D)] = 0



