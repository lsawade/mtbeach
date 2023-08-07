from .beach import Beachball
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import Axes3D


def plot_beachball(L):

    bb = Beachball(L, res=10, nosubdivide=False)

    fig = plt.figure(figsize=(3, 3))
    ax = Axes3D(fig, auto_add_to_figure=False,
                facecolor=None, computed_zorder=True)
    fig.add_axes(ax)
    ax.add_collection3d(Poly3DCollection(
        bb.vertices, edgecolors='k', linewidth=1.0,
        facecolors=np.where(bb.colors == 1, 'red', 'white'),
        alpha=0.5))
    ax.set_axis_off()
    limits = [-0.6, 0.6]
    ax.axes.set_xlim3d(*limits)
    ax.axes.set_ylim3d(*limits)
    ax.axes.set_zlim3d(*limits)
    ax.set_box_aspect([1.0, 1.0, 1.0])
    ax.view_init(elev=20., azim=0)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    return fig, ax
