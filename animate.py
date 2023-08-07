# %%

from mtbeach.plot import plot_beachball
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.collections as mc
from matplotlib import animation
# %%
# bb = Beachball([1, 19, -20], res=10, nosubdivide=False)


def animate_beachball(outfile='beachball.gif'):
    """Animate a beachball by rotating it around the z-axis at constant
    elevation and rate."""

    fig, ax = plot_beachball([1, 19, -20])

    def init():
        return fig,

    def animate(i):
        print(i)
        ax.view_init(elev=20., azim=2*i)
        # fig.savefig(f'frames/frame{i:>010d}.png', transparent=True, dpi=300)
        return fig,

    ani = animation.FuncAnimation(
        fig, animate, 90, init_func=init, interval=1/30, blit=True)

    ani.save(outfile, writer='imagemagick', dpi=70, fps=30,
             savefig_kwargs=dict(transparent=True))

# %%


animate_beachball()
