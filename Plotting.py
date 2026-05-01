import numpy as np
import matplotlib.pyplot as plt


X_p = np.loadtxt("X_p.csv",float,delimiter=',')
Y_p = np.loadtxt("Y_p.csv",float,delimiter=',')
X_u = np.loadtxt("X_u.csv",float,delimiter=',')
Y_u = np.loadtxt("Y_u.csv",float,delimiter=',')
X_v = np.loadtxt("X_v.csv",float,delimiter=',')
Y_v = np.loadtxt("Y_v.csv",float,delimiter=',')


t = 9.55

u = np.loadtxt(f"u_Time={t}.csv",float,delimiter=',')
v = np.loadtxt(f"v_Time={t}.csv",float,delimiter=',')
p = np.loadtxt(f"p_Time={t}.csv",float,delimiter=',')


plt.contourf(X_p,Y_p,p)

x_uf = X_u.ravel()
y_uf = Y_u.ravel()
u_f = u.ravel()

x_vf = X_v.ravel()
y_vf = Y_v.ravel()
v_f = v.ravel()


plt.quiver(x_uf, y_uf, u_f, np.zeros_like(u_f), scale=50)

plt.show()

plt.contourf(X_p,Y_p,p)

plt.quiver(x_vf, y_vf, np.zeros_like(v_f), v_f, scale=50)




plt.show()
