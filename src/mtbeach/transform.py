import numpy as np


def sphere2cart(r, theta, phi):
    """Mathematical transformation from spherical to cartesian coordinates.

    Parameters
    ----------
    r : arraylike
        Radial distance from the origin.
    phi : arraylike
        Angle from the z-axis.
    theta : arraylike
        Angle from the x-axis.

    Returns
    -------
    x : arraylike
        Cartesian x-coordinate.
    y : arraylike
        Cartesian y-coordinate.
    z : arraylike
        Cartesian z-coordinate.

    Notes
    -----
    r, theta, and phi must be the same shape.
    """
    x = r * np.sin(phi) * np.cos(theta)
    y = r * np.sin(phi) * np.sin(theta)
    z = r * np.cos(phi)

    return x, y, z


def cart2sphere(x, y, z):
    """Mathematical transformation from cartesian to spherical coordinates.

    Parameters
    ----------
    x : arraylike
        Cartesian x-coordinate.
    y : arraylike
        Cartesian y-coordinate.
    z : arraylike
        Cartesian z-coordinate.

    Returns
    -------
    r : arraylike
        Radial distance from the origin.
    theta : arraylike
        Angle from the z-axis.
    phi : arraylike
        Angle from the x-axis.

    Notes
    -----
    x, y, and z must be the same shape.
    """
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(z/r)
    phi = np.arctan2(y, x)

    return r, theta, phi


def cart2polar(x, y):
    """Mathematical transformation from cartesian to polar coordinates.

    Parameters
    ----------
    x : arraylike
        Cartesian x-coordinate.
    y : arraylike
        Cartesian y-coordinate.

    Returns
    -------
    r : arraylike
        Radial distance from the origin.
    theta : arraylike
        Angle from the z-axis.

    Notes
    -----
    x and y must be the same shape.
    """
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)

    return r, theta


def polar2cart(r, theta):
    """Mathematical transformation from polar to cartesian coordinates.

    Parameters
    ----------
    r : arraylike
        Radial distance from the origin.
    theta : arraylike
        Angle from the z-axis.

    Returns
    -------
    x : arraylike
        Cartesian x-coordinate.
    y : arraylike
        Cartesian y-coordinate.

    Notes
    -----
    r and theta must be the same shape.
    """
    x = r * np.cos(theta)
    y = r * np.sin(theta)

    return x, y
