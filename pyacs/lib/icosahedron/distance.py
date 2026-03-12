"""Euclidean distance between two 3D points."""

import math


def distance(A, B):
    """Return Euclidean distance between two points A and B (each [x, y, z])."""
    return math.sqrt((A[0] - B[0])**2 + (A[1] - B[1])**2 + (A[2] - B[2])**2)
