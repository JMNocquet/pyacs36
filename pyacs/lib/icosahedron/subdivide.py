"""Subdivide icosahedron faces and push vertices to unit sphere."""

from pyacs.lib.euclid import Vector3


def subdivide(verts, faces):
    """Subdivide each triangle into four triangles, pushing verts to the unit sphere.

    Parameters
    ----------
    verts : list
        List of vertices (modified in place).
    faces : list
        List of face index tuples (modified in place).

    Returns
    -------
    verts : list
        Updated vertices.
    faces : list
        Updated faces.
    """
    from progress.bar import Bar

    triangles = len(faces)
    bar = Bar('', max=triangles, suffix='%(percent).1f%% - %(elapsed)ds')
    for faceIndex in range(triangles):
        face = faces[faceIndex]
        a, b, c = (Vector3(*verts[vertIndex]) for vertIndex in face)
        verts.append((a + b).normalized()[:])
        verts.append((b + c).normalized()[:])
        verts.append((a + c).normalized()[:])

        i = len(verts) - 3
        j, k = i + 1, i + 2
        faces.append((i, j, k))
        faces.append((face[0], i, k))
        faces.append((i, face[1], j))
        faces[faceIndex] = (k, j, face[2])

        bar.next()

    bar.finish()

    return verts, faces
