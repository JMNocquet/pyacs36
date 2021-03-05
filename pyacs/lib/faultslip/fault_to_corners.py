
def fault_to_corners( lon, lat, depth, length , width, strike, dip):
    """

    :param lon: origin longitude for fault
    :param lat: origin latitude for fault
    :param length: fault length in km
    :param width:  fault width in km
    :param strike: fault strike in degrees
    :param dip: fault dip in degrees
    :return: lonc, latc, depthc 1D numpy arrays of the 4 fault corners

    :note: use a flat Earth approximation, depth is negative
    """

    # import
    import numpy as np
    import pyacs.lib.coordinates as coo
    import pyacs.message
    # convert to xy
    xref,yref = coo.geo2flat_earth(lon,lat)
    if depth >0:
        pyacs.message.warning("depth>0. Changed to negative.")
        depth = -depth


    # computes delta_x, delta_y, delta_z for the rectangular dislocation vertices
    alpha = np.radians(-strike) + np.pi / 2.

    delta_x = np.cos(np.radians( dip )) * np.cos(alpha - np.pi / 2.) * width
    delta_y = np.cos(np.radians( dip )) * np.sin(alpha - np.pi / 2.) * width
    delta_z = np.sin(np.radians( dip )) * width
    DELTA_BOTTOM = np.array([delta_x, delta_y, -delta_z]).T

    delta_x =  np.cos(alpha) * length
    delta_y =  np.sin(alpha) * length
    delta_z = delta_x * 0.0
    DELTA_TOP = np.array([delta_x, delta_y, delta_z]).T


    # get rectangle coordinates
    A = np.array([xref,yref,depth])
    B = A + DELTA_TOP
    C = B + DELTA_BOTTOM
    D = A + DELTA_BOTTOM

    # convert back to geographical coordinates
    A[0],A[1] = coo.flat_earth2geo(A[0],A[1])
    B[0],B[1] = coo.flat_earth2geo(B[0],B[1])
    C[0],C[1] = coo.flat_earth2geo(C[0],C[1])
    D[0],D[1] = coo.flat_earth2geo(D[0],D[1])

    # return
    return np.array([A,B,C,D])




