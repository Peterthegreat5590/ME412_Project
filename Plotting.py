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


p_mask = np.zeros_like(p, dtype=bool)
p_mask[(X_p<0)|(X_p>5)|(Y_p<0)|(Y_p>6)]=True
p_mask[(X_p>1)&(X_p<4)&(Y_p>1)] = True

p_masked = np.ma.array(p,mask=p_mask)
X_p_masked = np.ma.array(X_p, mask=p_mask)
Y_p_masked = np.ma.array(Y_p, mask=p_mask)

plt.contourf(X_p_masked,Y_p_masked,p_masked,cmap = "jet", corner_mask = False)

x_uf = X_p.ravel()
y_uf = Y_p.ravel()
u_f = u.ravel()

u_f = u_f[(x_uf>0)&(x_uf<5)&(y_uf>0)&(y_uf<6)]
x_uf, y_uf = x_uf[(x_uf>0)&(x_uf<5)&(y_uf>0)&(y_uf<6)], y_uf[(x_uf>0)&(x_uf<5)&(y_uf>0)&(y_uf<6)]

u_f = u_f[(x_uf<1)|(x_uf>4)|(y_uf<1)]
x_uf, y_uf = x_uf[(x_uf<1)|(x_uf>4)|(y_uf<1)], y_uf[(x_uf<1)|(x_uf>4)|(y_uf<1)]


x_vf = X_p.ravel()
y_vf = Y_p.ravel()
v_f = v.ravel()

v_f = v_f[(x_vf>0)&(x_vf<5)&(y_vf>0)&(y_vf<6)]
x_vf, y_vf = x_vf[(x_vf>0)&(x_vf<5)&(y_vf>0)&(y_vf<6)], y_vf[(x_vf>0)&(x_vf<5)&(y_vf>0)&(y_vf<6)]

v_f = v_f[(x_vf<1)|(x_vf>4)|(y_vf<1)]
x_vf, y_vf = x_vf[(x_vf<1)|(x_vf>4)|(y_vf<1)], y_vf[(x_vf<1)|(x_vf>4)|(y_vf<1)]

plt.quiver(x_uf, y_uf, u_f, v_f, scale=6)

plt.xlim([0,5])
plt.ylim([0,6])

ax = plt.gca()

ax.set_aspect('equal', adjustable='box')

dx = X_p[0,1] - X_p[0,0]
dy = Y_p[1,0] - Y_p[0,0]
p_in = p[(X_p==(1/2)+dx/2)&(Y_p==(6-dy/2))]
p_out = p[(X_p==(9/2)+dx/2)&(Y_p==(6-dy/2))]
print(p_in-p_out)

velo = np.sqrt(u**2+v**2)
print(np.max(velo))

plt.show()


# plt.contourf(X_p_masked,Y_p_masked,p_masked,cmap = "jet", corner_mask = False)

# plt.quiver(x_vf, y_vf, np.zeros_like(v_f), v_f, scale=100)

# plt.xlim([0,5])
# plt.ylim([0,6])


# plt.show()
