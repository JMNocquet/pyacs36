def make_distance_matrix_from_sgts(self):
    """Compute XYZ distance matrix (km) between all sites in this Sgts.

    Returns
    -------
    numpy.ndarray
        Distance matrix in km (same order as sorted lcode()).
    """
    import numpy as np
    import pyacs.lib.coordinates

    sgts = self.copy()

    def sgts_2_coo(sgts):
        coo = np.zeros((sgts.n(), 3))
        for i, code in enumerate(sorted(sgts.lcode())):
            coo[i, 0] = sgts.__dict__[code].lon
            coo[i, 1] = sgts.__dict__[code].lat
            coo[i, 2] = sgts.__dict__[code].h
        return coo

    def make_distance_matrix_fast(coo):
        coo_xyz = np.array(pyacs.lib.coordinates.geo2xyz(coo[:, 0], coo[:, 1], coo[:, 2], unit='dec_deg')).T
        x1x2 = np.matmul(coo_xyz[:, 0].reshape(-1, 1), coo_xyz[:, 0].reshape(-1, 1).T).T
        y1y2 = np.matmul(coo_xyz[:, 1].reshape(-1, 1), coo_xyz[:, 1].reshape(-1, 1).T).T
        z1z2 = np.matmul(coo_xyz[:, 2].reshape(-1, 1), coo_xyz[:, 2].reshape(-1, 1).T).T
        X1 = coo_xyz[:, 0] ** 2 + (coo_xyz[:, 0] ** 2).reshape(-1, 1) - 2 * x1x2
        Y1 = coo_xyz[:, 1] ** 2 + (coo_xyz[:, 1] ** 2).reshape(-1, 1) - 2 * y1y2
        Z1 = coo_xyz[:, 2] ** 2 + (coo_xyz[:, 2] ** 2).reshape(-1, 1) - 2 * z1z2
        return np.sqrt(X1 + Y1 + Z1) / 1000.

    coor = sgts_2_coo(sgts)
    return make_distance_matrix_fast(coor)
