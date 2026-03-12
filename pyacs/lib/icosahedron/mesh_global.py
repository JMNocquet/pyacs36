"""Build global icosahedron mesh by subdivision."""

from .build_icosahedron import icosahedron
from .subdivide import subdivide
from .distance import distance


def mesh_global(num_subdivisions=6):
    """Build global icosahedron mesh (verts and faces on unit sphere).

    Parameters
    ----------
    num_subdivisions : int, optional
        Number of subdivision iterations. Default is 6.

    Returns
    -------
    verts : list
        Vertices on unit sphere.
    faces : list
        Face index tuples.
    """
    Rt = 6371.0E3
    (verts, faces) = icosahedron()
    print("-- Number of vertices and faces of initial icosahedron ", len(verts), len(faces))
    print("-- Number of divisions ", num_subdivisions)
    print("-- Now doing subdivision...")
    for x in range(num_subdivisions):
        print("   - Division iteration: ", x + 1, "/", num_subdivisions)
        verts, faces = subdivide(verts, faces)
        print("   - New Number of vertices and faces ", len(verts), len(faces))
        A = verts[faces[0][0]]
        B = verts[faces[0][1]]
        print(("   - New triangle edge distance : %8.3lf km" % (distance(A, B) * Rt / 1000.)))
    return (verts, faces)
