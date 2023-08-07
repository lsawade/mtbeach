# %%
from mtbeach.poly import Triangle, Quadrilateral, Pentagon
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import PolyCollection
import typing as tp

# %%


def plot_split(poly: Quadrilateral, f: tp.Callable):
    # Add axes
    ax = plt.gca()

    # Get the colors
    colors = ['yellow', 'orange']

    # Split the polygon
    polys = poly.subdivide(f)

    # Plot the polygon
    ax.add_collection(PolyCollection([poly.vertices], edgecolors='k', linewidth=0.25,
                                     facecolors='none',))

    # Plot the split polygon
    ax.add_collection(PolyCollection([_poly.vertices for _poly in polys], edgecolors='k', linewidth=0.25,
                                     facecolors=colors[: len(polys)], alpha=0.75))

    for _poly in polys:
        print("polytype", type(_poly))
        if isinstance(_poly, Pentagon):
            print('Pentagon')
            xc, yc = _poly.COM
        elif isinstance(_poly, Quadrilateral):
            print('Rectangle')
            xc, yc = _poly.centroid
        elif isinstance(_poly, Triangle):
            print('Triangle')
            xc, yc = _poly.centroid
        else:
            continue

        print('COM(x,y):', xc, yc)
        ax.plot(xc, yc, 'kx', lw=2.0)

    x = np.linspace(-0.5, 2.5, 5)
    y = f(x)
    ax.plot(x, y, 'k', lw=2.0)

    x1, y1 = poly.vertices[0]
    x2, y2 = poly.vertices[1]
    x3, y3 = poly.vertices[2]
    x4, y4 = poly.vertices[3]

    kwargs = dict(ha='center', va='center', backgroundcolor='w')
    ax.text(x1-0.2, y1-0.2, '(x1,y1)', **kwargs)
    ax.text(x2+0.2, y2-0.2, '(x2,y2)', **kwargs)
    ax.text(x3+0.2, y3+0.2, '(x3,y3)', **kwargs)
    ax.text(x4-0.2, y4+0.2, '(x4,y4)', **kwargs)

    ax.set_xlim(-0.5, 2.5)
    ax.set_ylim(-0.25, 1.25)

    ax.axis('off')


x1, x4 = 0, 0
x2, x3 = 2, 2
y1, y2 = 0, 0
y3, y4 = 1, 1

point1 = [x1, y1]
point2 = [x2, y2]
point3 = [x3, y3]
point4 = [x4, y4]
vertices = np.array([point1, point2, point3, point4])


def f1(x): return 2 * x + 0.5
def f2(x): return 0.5 * x - 0.5
def f3(x): return 2.0 * x + 1.5


poly = Quadrilateral(vertices)


plt.figure(figsize=(9, 1.5))
plt.subplot(1, 3, 1)
plot_split(poly, f1)
plt.subplot(1, 3, 2)
plot_split(poly, f2)
plt.subplot(1, 3, 3)
plot_split(poly, f3)
plt.show()
