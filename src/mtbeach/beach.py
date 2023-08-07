import numpy as np
from .poly import BasicQuads, Polygons
from functools import partial


def HueBeach(v, M, ColorIn=-1, ColorOut=1):
    if np.dot(np.dot(M, v), v) > 0:
        return ColorOut
    else:
        return ColorIn


zref = np.diag([1, 1, -1])


def phiN(theta, L):
    """phiN is the phi coordinate of the (upper) nodal point at theta. phiN
    comes from solving (lambda1 x, lambda2, y, lambda3 z) .dot. (x,y,z) = 0,
    which is the nodal curve, which is good for any (lambda1, lambda2, lambda3).
    But to get phi for any theta you need lambda1 > 0, lambda2 > 0, lambda3 < 0
    OR lambda1 < 0, lambda2 < 0, lambda3 > 0.
    So, (theta, lambda1, lambda2, lambda3) -> (theta, phiN(theta))

    Parameters
    ----------
    theta : float or arraylike
        spherical longitude
    L : iterable
        lambda 1, 2, 3

    Returns
    -------
    _type_
        _description_
    """
    numerator = L[0] * np.cos(theta) ** 2 + \
        L[1] * np.sin(theta) ** 2
    denominator = L[0] * np.cos(theta) ** 2 + \
        L[1] * np.sin(theta) ** 2 - L[2]
    return np.arccos(np.sqrt(numerator/denominator))


class Beachball(Polygons):

    def __init__(self, L, res=10, nosubdivide=False):

        self.L = L
        self.res = res
        self.nosubdivide = nosubdivide

        self._calc_beachball()

        super().__init__(self.polys)

        del self.polys

        self.colors = self._calc_colors()

    def _calc_beachball(self):

        self.polys = BasicQuads(
            np.linspace(0, 2*np.pi, int((360/self.res)) + 1, endpoint=True),
            np.linspace(0, np.pi, int((180/self.res)) + 1, endpoint=True)
        )

        # For debugging the polygon division
        # self.polys = BasicQuads(
        #     np.linspace(0, 20/180*np.pi, 3, endpoint=True),
        #     np.linspace(70/180*np.pi, 90/180*np.pi, 3, endpoint=True))

        if self.nosubdivide is False:
            self._subdivide()

        # Make the quads spherical by adding a radial component
        self.polys.add_dim(1, dim=0)

        # Transform to cartesian coordinates
        self.polys.sphere2cart()

    def _calc_colors(self):
        """1 for compressional, -1 for dilatational, 0 for neutral"""
        colors = np.zeros(len(self))

        DL = np.diag(self.L)

        for _i, poly in enumerate(self):
            x = np.mean(poly.vertices, axis=0)
            colors[_i] = HueBeach(x, DL)

        return colors

    def _subdivide(self):

        if (((self.L[0] >= 0) and (self.L[1] >= 0) and (self.L[2] >= 0)) or
                ((self.L[0] <= 0) and (self.L[1] <= 0) and (self.L[2] <= 0))):

            return

        else:

            print('we are here')
            pphiN = partial(phiN, L=self.L)

            def tpphiN(theta):
                return np.pi - pphiN(theta)

            polys = []
            for _poly in self.polys:

                npoly = _poly.subdivide(pphiN)
                if len(npoly) == 2:
                    polys.extend(npoly)
                else:
                    polys.extend(_poly.subdivide(tpphiN))

            self.polys = Polygons(polys)
