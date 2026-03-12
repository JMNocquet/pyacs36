"""Construct an icosahedron on the unit sphere."""

import math


def icosahedron():
    """Construct an icosahedron on the unit sphere.

    Returns
    -------
    verts : list of tuple
        Vertices (x, y, z) on unit sphere.
    faces : list of tuple
        Face indices, each a 3-tuple of vertex indices.
    """
    tau = (1 + math.sqrt(5)) / 2
    radius = math.sqrt(1 + tau**2)

    verts = [
        (0, 1, tau),
        (0, 1, -tau),
        (0, -1, tau),
        (0, -1, -tau),
        (1, tau, 0),
        (1, -tau, 0),
        (-1, tau, 0),
        (-1, -tau, 0),
        (tau, 0, 1),
        (tau, 0, -1),
        (-tau, 0, 1),
        (-tau, 0, -1),
    ]

    faces = [
        (0, 2, 8),
        (0, 2, 10),
        (0, 4, 6),
        (0, 4, 8),
        (0, 6, 10),
        (1, 3, 9),
        (1, 3, 11),
        (1, 4, 6),
        (1, 4, 9),
        (1, 6, 11),
        (2, 5, 7),
        (2, 5, 8),
        (2, 7, 10),
        (3, 5, 7),
        (3, 5, 9),
        (3, 7, 11),
        (4, 8, 9),
        (5, 8, 9),
        (6, 10, 11),
        (7, 10, 11),
    ]

    verts = [
        (0, 1 / radius, tau / radius),
        (0, 1 / radius, -tau / radius),
        (0, -1 / radius, tau / radius),
        (0, -1 / radius, -tau / radius),
        (1 / radius, tau / radius, 0),
        (1 / radius, -tau / radius, 0),
        (-1 / radius, tau / radius, 0),
        (-1 / radius, -tau / radius, 0),
        (tau / radius, 0, 1 / radius),
        (tau / radius, 0, -1 / radius),
        (-tau / radius, 0, 1 / radius),
        (-tau / radius, 0, -1 / radius),
    ]
    return (verts, faces)
