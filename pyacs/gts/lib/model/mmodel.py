######################################################################################################################
## Generates a model time series
######################################################################################################################

def mmodel(self):
    """
    Generates a modeled time series from the parameters read in self
    """

    import numpy as np
    from pyacs.gts.Gts import Gts
    import inspect

    if isinstance(self.offsets_values, np.ndarray):
        noffset = len(self.offsets_values)
        offsets_values = np.zeros((noffset, 7))
    else:
        noffset = 0
        offsets_values = None

    t_ref = self.data[0, 0]
    t_ref_seasonal = 2010.0

    ### vel = [vel_N vel_E vel_U svel_N s_vel_E svel_U]
    vel = np.zeros(6)
    ### annual = [amplitude_NEU phase_NEU]
    annual = []
    ### semi_annual = [amplitude_NEU phase_NEU]
    semi_annual = []

    data = np.copy(self.data)

    residuals = np.zeros(data.shape)
    residuals[:, 0] = data[:, 0]

    ndate = len(data[:, 0])

    for k in range(1, 4):
        ## write matrix A in general case
        A = np.zeros([ndate, (6 + noffset)], float)
        for i in range(ndate):
            ti = data[i, 0]
            A[i, 0], A[i, 1], A[i, 2], A[i, 3], A[i, 4], A[i, 5] = 1., (ti - t_ref), \
                                                                   np.cos(2. * np.pi * (ti - t_ref_seasonal)), np.sin(
                2. * np.pi * (ti - t_ref_seasonal)), \
                                                                   np.cos(4. * np.pi * (ti - t_ref_seasonal)), np.sin(
                4. * np.pi * (ti - t_ref_seasonal))
            ## for offsets
            for j in range(noffset):
                if ti > self.offsets_dates[j]: A[i, (6 + j)] = 1.

        X = np.zeros([(6 + noffset)], float)

        v = self.velocity[k - 1]
        if isinstance(self.annual, np.ndarray):
            annual_p = self.annual[2 * (k - 1)]
            annual_q = self.annual[2 * (k - 1) + 1]
        else:
            annual_p = 0
            annual_q = 0

        if isinstance(self.semi_annual, np.ndarray):
            semi_annual_p = self.semi_annual[2 * (k - 1)]
            semi_annual_q = self.semi_annual[2 * (k - 1) + 1]
        else:
            semi_annual_p = 0
            semi_annual_q = 0

        X[0] = 0
        X[1] = v
        X[2] = annual_p
        X[3] = annual_q
        X[4] = semi_annual_p
        X[5] = semi_annual_q

        if isinstance(self.offsets_values, np.ndarray):
            for i in range(noffset):
                # print "self.offsets_values[i,k-1] ",self.offsets_values[i,k-1]
                X[6 + i] = self.offsets_values[i, k]

        else:
            for i in range(noffset):
                X[6 + i] = 0.0

        data[:, k] = np.dot(A, X)
    new_Gts = self.copy(data=data)

    return (new_Gts)
