# %%

from mtbeach.beach import Beachball, phiN
import numpy as np
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.collections as mc
import matplotlib.pyplot as plt
from mtbeach.poly import BasicQuads


theta = np.linspace(0, 2*np.pi, 37, endpoint=True)
phi = np.linspace(0, np.pi, 19, endpoint=True)

polys = BasicQuads(theta, phi)

# %%
plt.figure()
ax = plt.axes()
ax.add_collection(mc.PolyCollection(
    polys.vertices, edgecolor='k', facecolor='w', linewidth=0.5))
ax.set_xlim(0, 2*np.pi)
ax.set_ylim(0, np.pi)
plt.show()

# %%

# %%

# Make new set of polys
polys = BasicQuads(theta, phi)

# Make the quads spherical by adding a radial component
polys.add_dim(1, dim=0)

# Transform to cartesian coordinates
polys.sphere2cart()

fig = plt.figure()
ax = Axes3D(fig, facecolor=None, computed_zorder=True)
ax.add_collection3d(Poly3DCollection(
    polys.vertices, edgecolors='k', linewidth=0.1, facecolors='w', alpha=1))
ax.set_xlim(-1.1, 1.1)
ax.set_ylim(-1.1, 1.1)
ax.set_zlim(-1.1, 1.1)
ax.set_box_aspect([1.0, 1.0, 1.0])
plt.show()

# %%


def sphere2cartline(theta, phi):
    """Mathematical convention!!!"""
    xyz = np.zeros((len(theta), 3))
    xyz[:, 0] = np.sin(phi) * np.cos(theta)
    xyz[:, 1] = np.sin(phi) * np.sin(theta)
    xyz[:, 2] = np.cos(phi)

    return xyz


bb = Beachball([1, 19, -20], res=10, nosubdivide=False)
# theta = np.linspace(0, 2*np.pi, 1000)
# L = np.array([1, 19, -20])
# tL = np.diag([1, 1, -1]) @ L
# phi = phiN(theta, L)
# Uxyz = sphere2cartline(theta, phi)
# Lxyz = sphere2cartline(theta, np.pi-phi)
fig = plt.figure()
ax = Axes3D(fig, auto_add_to_figure=False,
            facecolor=None, computed_zorder=False)
fig.add_axes(ax)
scale = 1.025
# ax.plot(
#     Uxyz[:, 0]*scale,
#     Uxyz[:, 1]*scale,
#     Uxyz[:, 2]*scale,
#     'r', linewidth=5.0)
# ax.plot(
#     Lxyz[:, 0]*scale,
#     Lxyz[:, 1]*scale,
#     Lxyz[:, 2]*scale,
#     'b', linewidth=5.0)

ax.add_collection3d(Poly3DCollection(
    bb.vertices, edgecolors='k', linewidth=1.0,
    facecolors=np.where(bb.colors == 1, 'red', 'white'), alpha=1))
ax.set_xlim(-1.1, 1.1)
ax.set_ylim(-1.1, 1.1)
ax.set_zlim(-1.1, 1.1)
ax.set_box_aspect([1.0, 1.0, 1.0])
plt.show(block=False)
