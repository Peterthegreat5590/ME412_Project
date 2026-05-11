import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


X_p = np.loadtxt(Path("CSV/X_p.csv"),float,delimiter=',')
Y_p = np.loadtxt(Path("CSV/Y_p.csv"),float,delimiter=',')
X_u = np.loadtxt(Path("CSV/X_u.csv"),float,delimiter=',')
Y_u = np.loadtxt(Path("CSV/Y_u.csv"),float,delimiter=',')
X_v = np.loadtxt(Path("CSV/X_v.csv"),float,delimiter=',')
Y_v = np.loadtxt(Path("CSV/Y_v.csv"),float,delimiter=',')


t = 20

u = np.loadtxt(Path(f"CSV/u_Time={t:.3f}.csv"),float,delimiter=',')
v = np.loadtxt(Path(f"CSV/v_Time={t:.3f}.csv"),float,delimiter=',')
p = np.loadtxt(Path(f"CSV/p_Time={t:.3f}.csv"),float,delimiter=',')


u_prof = u[(X_u==5/2)&(Y_u<1.1)]
y_u_prof = Y_u[(X_u==5/2)&(Y_u<1.1)]

plt.plot(y_u_prof,u_prof)

plt.title("U velocity profile at X=2.5")
plt.xlabel("y")
plt.ylabel("U")
plt.xlim([0,1])
plt.ylim([0,plt.ylim()[1]])

ax = plt.gca()

ax.set_aspect('equal', adjustable='box')


plt.show()

v_prof = v[(Y_v==3.0)&(X_v<1.1)]
x_v_prof = X_v[(Y_v==3.0)&(X_v<1.1)]

plt.plot(x_v_prof,v_prof)

plt.title("V velocity profile at Y=3, X<1")
plt.xlabel("x")
plt.ylabel("V")
plt.xlim([0,1])
plt.ylim([plt.ylim()[0],0])

ax = plt.gca()

ax.set_aspect('equal', adjustable='box')


plt.show()

v_prof = v[(Y_v==3.0)&(X_v>3.9)]
x_v_prof = X_v[(Y_v==3.0)&(X_v>3.9)]

plt.plot(x_v_prof,v_prof)

plt.title("V velocity profile at Y=3, X>4")
plt.xlabel("x")
plt.ylabel("V")
plt.xlim([4,5])
plt.ylim([0,plt.ylim()[1]])

ax = plt.gca()

ax.set_aspect('equal', adjustable='box')


plt.show()