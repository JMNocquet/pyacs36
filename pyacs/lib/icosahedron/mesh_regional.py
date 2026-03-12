"""Build regional icosahedron mesh within bounds."""

import math

from pyacs.lib import coordinates
from shapely.geometry import Polygon as shp_polygon

from .build_icosahedron import icosahedron
from .subdivide import subdivide
from .distance import distance


def mesh_regional(num_subdivisions=6, bounds=None):
    """Build an equilateral triangle mesh on the unit sphere within bounds.

    Uses successive subdivisions of an icosahedron; only faces intersecting
    the given bounds are kept.

    Parameters
    ----------
    num_subdivisions : int, optional
        Number of subdivision iterations. Default is 6.
    bounds : str, optional
        Region as string, e.g. '/lon_min/lon_max/lat_min/lat_max' (degrees).

    Returns
    -------
    verts : list
        Vertices as [x, y, z] in geocentric cartesian on unit sphere.
    faces : list
        Faces as lists of vertex indices; faces[i][j] is vertex j of face i.
    """
    (lon_min, lon_max, lat_min, lat_max) = list(map(float, bounds.split('/')[1:]))
    (rlon_min, rlon_max, rlat_min, rlat_max) = list(map(math.radians, (lon_min, lon_max, lat_min, lat_max)))

    shp_rectangle = shp_polygon([
        (rlon_min, rlat_min), (rlon_max, rlat_min),
        (rlon_max, rlat_max), (rlon_min, rlat_max),
    ])

    Rt = 6371.0E3
    (verts, faces) = icosahedron()
    print("-- Number of vertices and faces of initial icosahedron ", len(verts), len(faces))
    print("-- Number of divisions ", num_subdivisions)
    print("-- Now doing subdivision...")
    ndivision = 6
    if num_subdivisions < 6:
        ndivision = num_subdivisions
    for x in range(ndivision):
        print("   - Division iteration: ", x + 1, "/", num_subdivisions)
        verts, faces = subdivide(verts, faces)
        print("   - New Number of vertices and faces ", len(verts), len(faces))
        A = verts[faces[0][0]]
        B = verts[faces[0][1]]
        print(("   - New triangle edge distance : %8.3lf km" % (distance(A, B) * Rt / 1000.)))

    print("-- Now doing selection on bounds ", bounds)
    lfaces = []
    for face in faces:
        lgeo = []
        for j in [0, 1, 2]:
            (x, y, z) = verts[face[j]]
            (rlon, rlat, rh) = coordinates.xyz2geospheric(x, y, z)
            lgeo.append((rlon, rlat))

        shp_triangle = shp_polygon(lgeo)
        shp_intersection = shp_rectangle.intersects(shp_triangle)
        if not shp_intersection:
            continue
        else:
            lfaces.append(face)

    print("-- Keeping ", len(lfaces), " faces")
    faces = lfaces

    print("-- Now dealing with finer subdivisions")
    for x in range(ndivision, num_subdivisions):
        print("   - Division iteration: ", x + 1, "/", num_subdivisions)
        verts, faces = subdivide(verts, faces)
        print("   - New Number of vertices and faces ", len(verts), len(faces))
        A = verts[faces[0][0]]
        B = verts[faces[0][1]]
        print(("   - New triangle edge distance : %8.3lf km" % (distance(A, B) * Rt / 1000.)))
        print("   - Selecting faces on bounds ", bounds)

        lfaces = []
        for face in faces:
            lgeo = []
            for j in [0, 1, 2]:
                (x, y, z) = verts[face[j]]
                (rlon, rlat, rh) = coordinates.xyz2geospheric(x, y, z)
                lgeo.append((rlon, rlat))

            shp_triangle = shp_polygon(lgeo)
            shp_intersection = shp_rectangle.intersects(shp_triangle)

            if not shp_intersection:
                continue
            lfaces.append(face)

    faces = lfaces
    print("-- ", len(faces), " faces kept.")

    return (verts, faces)
