"""Polygon class and related functions.

Note that all of these functions assume convex polygons!"""

import numpy as np
import typing as tp
import numpy.typing as npt
from .transform import sphere2cart, cart2sphere, cart2polar, polar2cart


class Polygon(object):

    vertices: npt.ArrayLike
    N: int
    Ndim: int

    def __init__(self, vertices):
        """Vertices must be a 2D array with shape (N, M) where N is the number of vertices and M is the number of dimensions.

        Parameters
        ----------
        vertices : arraylike
            Coordinates of the vertices of the polygon.
        """
        self.vertices = vertices
        self.N = self.vertices.shape[0]
        self.Ndim = self.vertices.shape[1]

    def COM(self):
        pass

    def polar2cart(self):
        if self.Ndim != 2:
            raise ValueError(
                'For conversion to cartesian coordinates the polygon must be 2D.')
        else:
            self.vertices = np.vstack(
                polar2cart(self.vertices[:, 0],
                           self.vertices[:, 1])).T

    def sphere2cart(self):
        if self.Ndim != 3:
            raise ValueError(
                'For conversion to spherical coordinates the polygon must be 3D.')
        else:
            self.vertices = np.vstack(
                sphere2cart(self.vertices[:, 0],
                            self.vertices[:, 1],
                            self.vertices[:, 2])).T

    def cart2sphere(self):
        if self.Ndim != 3:
            raise ValueError(
                'For conversion to spherical coordinates the polygon must be 3D.')
        else:
            self.vertices = np.vstack(
                cart2sphere(self.vertices[:, 0],
                            self.vertices[:, 1],
                            self.vertices[:, 2])).T

    def cart2polar(self):
        if self.Ndim != 2:
            raise ValueError(
                'For conversion to polar coordinates the polygon must be 2D.')
        else:
            self.vertices = np.vstack(
                cart2polar(self.vertices[:, 0],
                           self.vertices[:, 1])).T

    def add_dim(self, v: npt.ArrayLike, dim: int = 0):
        """Convert 2D array to 3D array by adding a column of zeros. Uses
        numpys insert function dim can be 0 to Ndim where zero prepends a column
        and dim=Ndim appends a column.
        """

        if (v.shape[0] != self.N):
            raise ValueError(
                'Must have same number of elements as vertices in the polygon.')

        if (len(v.shape) != 1):
            raise ValueError('Must be 1D array.')

        # Inserts values
        self.vertices = np.insert(self.vertices, dim, v, axis=1)

        self.Ndim += 1

    def remove_dim(self, dim: int):
        """Dimension reduction from 3D to 2D by removing a column."""

        self.vertices = np.delete(self.vertices, dim, axis=1)
        self.Ndim -= 1


class Triangle(Polygon):

    def __init__(self, vertices):
        super().__init__(vertices)
        if self.N != 3:
            raise ValueError('Number of vertices must be 3.')

    @property
    def area(self):
        """Heron's formula"""

        a = np.linalg.norm(self.vertices[1, :] - self.vertices[0, :])
        b = np.linalg.norm(self.vertices[2, :] - self.vertices[1, :])
        c = np.linalg.norm(self.vertices[0, :] - self.vertices[2, :])

        s = 0.5 * (a + b + c)

        return np.sqrt(s * (s - a) * (s - b) * (s - c))

    @property
    def centroid(self):
        """Centroid of a triangle"""
        return np.mean(self.vertices, axis=0)


class Quadrilateral(Polygon):

    def __init__(self, vertices):
        super().__init__(vertices)
        if self.N != 4:
            raise ValueError('Number of vertices must be 4.')

    @property
    def centroid(self):
        """Centroid of a quadrilateral"""
        s1 = np.array([0, 1, 2])
        s2 = np.array([2, 3, 0])

        T1 = Triangle(self.vertices[s1, :])
        T2 = Triangle(self.vertices[s2, :])

        c1 = T1.centroid
        c2 = T2.centroid

        A1 = T1.area
        A2 = T2.area

        return (A1*c1 + A2*c2) / (A1 + A2)

    def subdivide(self, f: tp.Callable):
        """Subdivide a quadrilateral into two triangles.

        Parameters
        ----------
        f : callable
            Function of x that intersects with quad.

        The quad is aligned like so
        .. code::

            4---3
            |   |
            1---2

        Parameters
        ----------
        p: Quadrilateral
        2 floats
        f: callable
        function of x that intersects with quad

        Returns
        -------
        Set of polygons

        """

        # Corner points of the quad
        x1, y1 = self.vertices[0, :]
        x2, y2 = self.vertices[1, :]
        x3, y3 = self.vertices[2, :]
        x4, y4 = self.vertices[3, :]

        debug = False

        if debug:
            print(f'x1, y1 = {x1, y1}')
            print(f'x2, y2 = {x2, y2}')
            print(f'x3, y3 = {x3, y3}')
            print(f'x4, y4 = {x4, y4}')
            print(f'f(x1), f(x2), f(x3), f(x4) = {f(x1), f(x2), f(x3), f(x4)}')

        # Cuts quad's left and right sides
        # 4---3
        # x   x
        # 1---2
        if (((y1 <= f(x1)) and (f(x1) <= y4)) and
                ((y2 <= f(x2)) and (f(x2) <= y3))):

            if debug:
                print("Cuts quad's left and right sides")
            # Lower quad
            poly1 = Quadrilateral(np.array([[x1, y1],
                                            [x2, y2],
                                            [x2, f(x2)],
                                            [x1, f(x1)]]))
            # Upper quad
            poly2 = Quadrilateral(np.array([[x4, y4],
                                            [x3, y3],
                                            [x2, f(x2)],
                                            [x1, f(x1)]]))

        # Cuts quad's left and bottom side
        # 4---3
        # x   |
        # 1-x-2
        elif (((y1 <= f(x1)) and (f(x1) <= y4)) and
              (f(x2) <= y2)):
            if debug:
                print("Cuts quad's left and bottom side")

            # Lower left triangle
            poly1 = Triangle(np.array(
                [[x1, f(x1)],
                 [x1, y1],
                 [x1 + ((y1 - f(x1)) * (x2 - x1))/(f(x2) - f(x1)), y1]]))

            poly2 = Pentagon(np.array(
                [[x1, f(x1)],
                 [x1 + ((y1 - f(x1)) * (x2 - x1))/(f(x2) - f(x1)), y1],
                 [x2, y2],
                 [x3, y3],
                 [x4, y4]]))

        # Cuts quad's bottom and right side
        # 4---3
        # |   x
        # 1-x-2
        elif ((f(x1) <= y1) and
              ((y2 <= f(x2)) and (f(x2) <= y3))):
            if debug:
                print("Cuts quad's bottom and right side")

            # Lower right triangle
            poly1 = Triangle(np.array(
                [[x1 + ((y1 - f(x1)) * (x2 - x1))/(f(x2) - f(x1)), y1],
                 [x2, y2],
                 [x2, f(x2)]]))

            # Remaining pentagon
            poly2 = Pentagon(np.array(
                [[x1 + ((y1 - f(x1)) * (x2 - x1))/(f(x2) - f(x1)), y1],
                 [x2, f(x2)],
                 [x3, y3],
                 [x4, y4],
                 [x1, y1]]))

        # Cuts top and right side of quad:
        # 4-x-3
        # |   x
        # 1---2
        elif ((y4 <= f(x1)) and
                ((y2 <= f(x2)) and (f(x2) <= y3))):
            if debug:
                print("Cuts top and right side of quad")

            # Upper right triangle
            poly1 = Triangle(np.array(
                [[x1 + ((y3 - f(x1)) * (x2 - x1))/(f(x2) - f(x1)), y3],
                 [x3, y3],
                 [x3, f(x3)]]))

            # Remaining pentagon
            poly2 = Pentagon(np.array(
                [[x1 + ((y3 - f(x1)) * (x2 - x1))/(f(x2) - f(x1)), y3],
                 [x3, f(x3)],
                 [x2, y2],
                 [x1, y1],
                 [x4, y4]]))

        # Cuts left and top side of quad:
        # 4-x-3
        # x   |
        # 1---2
        elif (((y1 <= f(x1)) and (f(x1) <= y4)) and
                (y3 <= f(x2))):
            if debug:
                print("Cuts left and top side of quad")

            # Upper left triangle
            poly1 = Triangle(np.array(
                [[x1 + ((y3 - f(x1)) * (x2 - x1))/(f(x2) - f(x1)), y3],
                 [x4, y4],
                 [x1, f(x1)]]))

            # Remaining pentagon
            poly2 = Pentagon(np.array(
                [[x1 + ((y3 - f(x1)) * (x2 - x1))/(f(x2) - f(x1)), y3],
                 [x3, y3],
                 [x2, y2],
                 [x1, y1],
                 [x1, f(x1)]]))

        # Cuts top and bottom side of quad:
        # 4-x-3
        # |   |
        # 1-x-2
        elif ((f(x2) <= y2) and (y4 <= f(x4)) or
                (f(x1) <= y1) and (y3 <= f(x3))):
            if debug:
                print("Cuts top and bottom side of quad")

            # Left quad
            poly1 = Quadrilateral(np.array(
                [[x4 + ((y1-f(x4))*(x2-x4))/(f(x2)-f(x4)), y1],
                 [x2, y2],
                 [x2, y3],
                 [x4 + ((y4-f(x4))*(x2-x4))/(f(x2)-f(x4)), y3]]))

            # Right quad
            poly2 = Quadrilateral(np.array(
                [[x4 + ((y1-f(x4))*(x2-x4))/(f(x2)-f(x4)), y1],
                 [x1, y1],
                 [x4, y4],
                 [x4 + ((y4-f(x4))*(x2-x4))/(f(x2)-f(x4)), y3]]))

        else:
            return [self,]

        return [poly1, poly2]


class Pentagon(Polygon):

    def __init__(self, vertices):
        super().__init__(vertices)
        if self.N != 5:
            raise ValueError('Number of vertices must be 5.')

    @property
    def COM(self):
        s1 = np.array([0, 1, 4])
        s2 = np.array([1, 3, 4])
        s3 = np.array([1, 2, 3])

        T1 = Triangle(self.vertices[s1, :])
        T2 = Triangle(self.vertices[s2, :])
        T3 = Triangle(self.vertices[s3, :])

        A1 = T1.area
        A2 = T2.area
        A3 = T3.area
        c1 = T1.centroid
        c2 = T2.centroid
        c3 = T3.centroid

        return (A1*c1 + A2*c2 + A3*c3) / (A1 + A2 + A3)


class Polygons(list):

    def __init__(self, polys: tp.List[Polygon]):
        super().__init__(polys)

    def sphere2cart(self):
        for poly in self:
            poly.sphere2cart()

    def cart2sphere(self):
        for poly in self:
            poly.cart2sphere()

    def cart2polar(self):
        for poly in self:
            poly.cart2polar()

    def polar2cart(self):
        for poly in self:
            poly.polar2cart()

    def add_dim(self, v: float, dim: int = 0):
        """Adds a dimension to all polygons in the list.
        with a constant value v."""
        for poly in self:
            poly.add_dim(v*np.ones(poly.N), dim)

    def remove_dim(self, dim: int):
        for poly in self:
            poly.remove_dim(dim)

    @property
    def vertices(self):
        vertlist = []
        for poly in self:
            vertlist.append(poly.vertices)
        return vertlist


def BasicQuads(rangex, rangey):
    NX = len(rangex)
    NY = len(rangey)
    polys = []
    for i in range(NX-1):
        for j in range(NY-1):
            polys.append(Quadrilateral(np.array(
                [[rangex[i],   rangey[j]],
                 [rangex[i+1], rangey[j]],
                 [rangex[i+1], rangey[j+1]],
                 [rangex[i],   rangey[j+1]]])))
    return Polygons(polys)
