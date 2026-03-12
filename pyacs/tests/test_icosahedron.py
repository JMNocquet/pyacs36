"""
Unit tests for pyacs.lib.icosahedron (mesh construction and subdivision).
"""

import math

import numpy as np
import pytest

import pyacs.lib.icosahedron as ICO


# -----------------------------------------------------------------------------
# distance
# -----------------------------------------------------------------------------


def test_distance_origin():
    """Euclidean distance from origin to point."""
    A = (0.0, 0.0, 0.0)
    B = (3.0, 4.0, 0.0)
    assert ICO.distance(A, B) == 5.0


def test_distance_same_point():
    """Distance from a point to itself is zero."""
    A = (1.0, -2.0, 3.0)
    assert ICO.distance(A, A) == 0.0


def test_distance_symmetric():
    """distance(A, B) == distance(B, A)."""
    A = (1.0, 0.0, 0.0)
    B = (0.0, 1.0, 0.0)
    assert ICO.distance(A, B) == ICO.distance(B, A)
    assert np.isclose(ICO.distance(A, B), math.sqrt(2.0), atol=1e-10)


# -----------------------------------------------------------------------------
# icosahedron (build_icosahedron)
# -----------------------------------------------------------------------------


def test_icosahedron_returns_verts_faces():
    """icosahedron() returns (verts, faces) with correct lengths."""
    verts, faces = ICO.icosahedron()
    assert len(verts) == 12
    assert len(faces) == 20


def test_icosahedron_vertices_on_unit_sphere():
    """All vertices lie on the unit sphere (norm = 1)."""
    verts, _ = ICO.icosahedron()
    for v in verts:
        x, y, z = v[0], v[1], v[2]
        norm = math.sqrt(x * x + y * y + z * z)
        assert np.isclose(norm, 1.0, atol=1e-10)


def test_icosahedron_faces_valid_indices():
    """Each face has three valid vertex indices."""
    verts, faces = ICO.icosahedron()
    n = len(verts)
    for face in faces:
        assert len(face) == 3
        for i in face:
            assert 0 <= i < n


def test_icosahedron_face_triangles_distinct():
    """Face vertex indices are distinct."""
    _, faces = ICO.icosahedron()
    for face in faces:
        assert len(set(face)) == 3


# -----------------------------------------------------------------------------
# subdivide (requires progress.bar)
# -----------------------------------------------------------------------------


def test_subdivide_one_iteration():
    """One subdivision: 20 faces -> 80 faces, 12 verts -> 72 verts."""
    pytest.importorskip("progress.bar")
    verts, faces = ICO.icosahedron()
    verts = list(verts)
    faces = list(faces)
    v2, f2 = ICO.subdivide(verts, faces)
    assert len(f2) == 80
    assert len(v2) == 72


def test_subdivide_vertices_on_unit_sphere():
    """After subdivision, all vertices remain on unit sphere."""
    pytest.importorskip("progress.bar")
    verts, faces = ICO.icosahedron()
    verts = list(verts)
    faces = list(faces)
    v2, _ = ICO.subdivide(verts, faces)
    for v in v2:
        x, y, z = v[0], v[1], v[2]
        norm = math.sqrt(x * x + y * y + z * z)
        assert np.isclose(norm, 1.0, atol=1e-10)


# -----------------------------------------------------------------------------
# mesh_global
# -----------------------------------------------------------------------------


def test_mesh_global_zero_subdivisions():
    """mesh_global(0) returns base icosahedron (12 verts, 20 faces)."""
    verts, faces = ICO.mesh_global(num_subdivisions=0)
    assert len(verts) == 12
    assert len(faces) == 20


def test_mesh_global_one_subdivision():
    """mesh_global(1) returns 72 verts and 80 faces."""
    pytest.importorskip("progress.bar")
    verts, faces = ICO.mesh_global(num_subdivisions=1)
    assert len(verts) == 72
    assert len(faces) == 80


def test_mesh_global_vertices_unit_sphere():
    """mesh_global(0) vertices are on unit sphere."""
    verts, _ = ICO.mesh_global(num_subdivisions=0)
    for v in verts:
        x, y, z = v[0], v[1], v[2]
        norm = math.sqrt(x * x + y * y + z * z)
        assert np.isclose(norm, 1.0, atol=1e-10)


# -----------------------------------------------------------------------------
# mesh_regional (requires shapely)
# -----------------------------------------------------------------------------


def test_mesh_regional_zero_subdivisions_global_bounds():
    """mesh_regional(0, global bounds) returns 12 verts and 20 faces."""
    pytest.importorskip("shapely")
    verts, faces = ICO.mesh_regional(num_subdivisions=0, bounds="/-180/180/-90/90")
    assert len(verts) == 12
    assert len(faces) == 20


def test_mesh_regional_vertices_unit_sphere():
    """mesh_regional vertices are on unit sphere."""
    pytest.importorskip("shapely")
    verts, _ = ICO.mesh_regional(num_subdivisions=0, bounds="/-180/180/-90/90")
    for v in verts:
        x, y, z = v[0], v[1], v[2]
        norm = math.sqrt(x * x + y * y + z * z)
        assert np.isclose(norm, 1.0, atol=1e-10)


def test_mesh_regional_restricted_bounds_fewer_faces():
    """Restricted bounds yields fewer faces than global."""
    pytest.importorskip("shapely")
    _, faces_global = ICO.mesh_regional(num_subdivisions=0, bounds="/-180/180/-90/90")
    _, faces_regional = ICO.mesh_regional(num_subdivisions=0, bounds="/0/30/40/60")
    assert len(faces_regional) < len(faces_global)
    assert len(faces_regional) >= 1
