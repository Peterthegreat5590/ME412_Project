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


mask = np.zeros_like(u, dtype=bool)
mask[(X_p<0)|(X_p>5)|(Y_p<0)|(Y_p>6)]=True
mask[(X_p>1)&(X_p<4)&(Y_p>1)] = True

V = np.sqrt(u**2+v**2)
V_masked = np.ma.array(V,mask=mask)
X_p_masked = np.ma.array(X_p, mask=mask)
Y_p_masked = np.ma.array(Y_p, mask=mask)


plt.contourf(X_p_masked,Y_p_masked,V_masked,cmap = "jet", corner_mask = False, levels = 50)

plt.title("Velocity Magnitude")

plt.xlim([0,5])
plt.ylim([0,6])

ax = plt.gca()

ax.set_aspect('equal', adjustable='box')


plt.show()
