import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata


X_p = np.loadtxt("X_p.csv", float, delimiter=',')
Y_p = np.loadtxt("Y_p.csv", float, delimiter=',')
X_u = np.loadtxt("X_u.csv", float, delimiter=',')
Y_u = np.loadtxt("Y_u.csv", float, delimiter=',')
X_v = np.loadtxt("X_v.csv", float, delimiter=',')
Y_v = np.loadtxt("Y_v.csv", float, delimiter=',')


t = 0.1

u = np.loadtxt(f"u_Time={t}.csv", float, delimiter=',')
v = np.loadtxt(f"v_Time={t}.csv", float, delimiter=',')
p = np.loadtxt(f"p_Time={t}.csv", float, delimiter=',')


fig, axes = plt.subplots(1, 1, figsize=(14, 6))




# -------------------------------------------------------
# New plot: smooth contourf + streamplot
# -------------------------------------------------------
plt.sca(axes)
axes.set_title("Smooth: tricontourf + streamplot")

# --- Smooth pressure contour using tricontourf on unmasked points ---
x_pf = X_p.ravel()
y_pf = Y_p.ravel()
p_f = p.ravel()

# Apply same domain + geometry mask as before (keep only valid fluid region)
fluid = (
    (x_pf > 0) & (x_pf < 5) & (y_pf > 0) & (y_pf < 6) &
    ~((x_pf > 1) & (x_pf < 4) & (y_pf > 1))
)
x_pf, y_pf, p_f = x_pf[fluid], y_pf[fluid], p_f[fluid]

# tricontourf interpolates scattered points smoothly without needing a regular grid
axes.tricontourf(x_pf, y_pf, p_f, levels=50, cmap="jet")

# --- Streamplot (curved velocity lines) on a regular grid ---
# Build a regular grid over the domain
nx, ny = 200, 240
xi = np.linspace(0, 5, nx)
yi = np.linspace(0, 6, ny)
Xi, Yi = np.meshgrid(xi, yi)

# Collect scattered u, v points (same fluid mask as quiver above)
x_vel = X_p.ravel()
y_vel = Y_p.ravel()
u_sc = u.ravel()
v_sc = v.ravel()

fluid_vel = (
    (x_vel > 0) & (x_vel < 5) & (y_vel > 0) & (y_vel < 6) &
    ~((x_vel > 1) & (x_vel < 4) & (y_vel > 1))
)
x_vel, y_vel = x_vel[fluid_vel], y_vel[fluid_vel]
u_sc, v_sc = u_sc[fluid_vel], v_sc[fluid_vel]

# Interpolate scattered velocity onto regular grid
Ui = griddata((x_vel, y_vel), u_sc, (Xi, Yi), method='linear')
Vi = griddata((x_vel, y_vel), v_sc, (Xi, Yi), method='linear')

# Mask the solid region (x: 1–4, y > 1) and outside domain on the grid
# Use a slightly expanded mask so griddata fill doesn't bleed into the block
solid = ((Xi >= 1) & (Xi <= 4) & (Yi >= 1)) | \
        (Xi <= 0) | (Xi >= 5) | (Yi <= 0) | (Yi >= 6)
Ui[solid] = np.nan
Vi[solid] = np.nan

speed = np.sqrt(Ui**2 + Vi**2)
axes.streamplot(Xi, Yi, Ui, Vi, color=speed, cmap="cool",
                   linewidth=1, density=1.5, arrowsize=1.2)

# Cover any streamline fragments that bled into the solid block with a patch
from matplotlib.patches import Rectangle
axes.add_patch(Rectangle((1, 1), 3, 5, color='white', zorder=5))

plt.xlim([0, 5])
plt.ylim([0, 6])
axes.set_aspect('equal', adjustable='box')

plt.tight_layout()
plt.show()