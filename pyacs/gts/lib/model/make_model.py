######################################################################################################################
## estimate the offsets of time serie: y = a + b.t + [c.sin(2.pi.t)+d.cos(2.pi.t)]+[e.sin(4.pi.t)+f.cos(4.pi.t)] + offsets
######################################################################################################################
def make_model(self, option='detrend', method='L2', loutlier=None, in_place=False):
    """
        Estimate linear model parameters using least squares
        input: data: Gts format
        option are: 'detrend'/'detrend_annual'/'detrend_seasonal'
        output: new Gts object: time series is now the residuals wrt to the model and its associated values (vel, annual, semi-annual etc)
    """

    import numpy as np
    from pyacs.gts.Gts import Gts
    import inspect

    ###########################################################################
    # check data is not None
    from pyacs.gts.lib.errors import GtsInputDataNone

    try:
        if self.data is None:
            # raise exception
            raise GtsInputDataNone(inspect.stack()[0][3], __name__, self)
    except GtsInputDataNone as error:
        # print PYACS WARNING
        print(error)
        return (self)
    ###########################################################################

    import pyacs.lib.glinalg as glinalg
    import pyacs.lib.robustestimators as RobEst

    #    from gts_estimators import least_square

    data = np.copy(self.data)

    ### REMOVE OUTLIERS FOR LINEAR ESTIMATION
    if loutlier: data = np.delete(data, loutlier, axis=0)

    ### REMOVE OFFSET DATES OUTSIDE THE TIME SPAN OF THE TIME SERIES

    if self.offsets_dates is not None and len(self.offsets_dates) > 0:
        # keep offsets_dates only within the time series
        sel_offsets_dates = \
            [self.offsets_dates[i] for i in range(len(self.offsets_dates)) if
             (self.offsets_dates[i] > self.data[0, 0] and self.offsets_dates[i] < self.data[-1, 0])]
        noffset = len(sel_offsets_dates)

        ### offsets_values = [time_offsets offsets_NEU s_offsets_NEU]
        offsets_values = np.zeros((noffset, 7))
        offsets_values[:, 0] = self.offsets_dates
    else:
        noffset = 0
        offsets_values = None

    # REF DATES

    t_ref = self.data[0, 0]
    t_ref_seasonal = 2010.0

    # INDEX
    if option == 'detrend_seasonal':
        del_index = []
        n_default_unknown = 6
    elif option == 'detrend_annual':
        del_index = [4, 5]
        n_default_unknown = 4
    elif option == 'detrend':
        del_index = [2, 3, 4, 5]
        n_default_unknown = 2
    else:
        print('    ERROR!!! check the option of estimation: detrend/detrend_seasonal/detrend_annual')

    # INIT VEL, ANNUAL, SEMI_ANNUAL, RESIDUALS ARRAYS
    ### vel = [vel_N vel_E vel_U svel_N s_vel_E svel_U]
    vel = np.zeros(6)
    ### annual = [amplitude_NEU phase_NEU]
    annual = []
    ### semi_annual = [amplitude_NEU phase_NEU]
    semi_annual = []

    residuals = np.zeros(data.shape)
    residuals[:, 0] = data[:, 0]

    ndate = len(data[:, 0])

    # BUILD LINEAR SYSTEM

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

        ### take the design matrix
        A = np.delete(A, del_index, axis=1)

        # solve

        (X, COV, V) = glinalg.lsw_full(A, data[:, k], data[:, k + 3])

        if method == 'L1':
            (X, V) = RobEst.Dikin(A, data[:, k], data[:, k + 3], eps=1.0E-4)

        s_X = np.sqrt(np.diag(COV))

        ## calculate residuals, predicted mb_files
        residuals[:, k] = V
        residuals[:, k + 3] = data[:, k + 3]

        ## velocity
        vel[k - 1] = X[1]
        vel[k + 2] = s_X[1]

        ## calculate the offset amplitudes
        if noffset > 0:
            offsets_values[:, k] = X[n_default_unknown:(n_default_unknown + noffset)]
            offsets_values[:, k + 3] = s_X[n_default_unknown:(n_default_unknown + noffset)]

        if option == 'detrend_annual':
            ## calculate the annual motion: amplitude & phase
            annual.append((X[2], X[3]))

        if option == 'detrend_seasonal':
            ## calculate the annual motion: amplitude & phase
            annual.append((X[2], X[3]))
            ## calculate the annual motion: amplitude & phase
            semi_annual.append((X[4], X[5]))

    if len(annual) > 0:
        annual = np.hstack(annual)
    else:
        annual = None
    if len(semi_annual) > 0:
        semi_annual = np.hstack(semi_annual)
    else:
        semi_annual = None

    new_Gts = self.copy()

    new_Gts.offsets_values = offsets_values
    new_Gts.annual = annual
    new_Gts.semi_annual = semi_annual
    new_Gts.velocity = vel

    new_Gts.data = residuals

    if in_place:
        self = new_Gts.copy()
        del new_Gts
        return (self)
    else:
        return (new_Gts)
