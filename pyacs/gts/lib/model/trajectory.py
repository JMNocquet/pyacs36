"""
Non linear trajectory models for Geodetic Time Series
"""


###############################################################################
def trajectory(self,
               model_type,
               offset_dates=[],
               eq_dates=[],
               H_fix={},
               H_constraints={},
               H_bounds={},
               component='NEU',
               verbose=False):
###############################################################################
    """
    Calculates the parameters of a (non-linear) trajectory model for a Geodetic Time Series.
    The trajectory model is:

    y(t) =

    trend : trend_cst + trend * ( t - t0 ) +

    annual: a_annual * cos( 2*pi + phi_annual ) +

    semi-annual: a_semi_annual * cos( 2*pi + phi_semi_annual ) +

    offset : Heaviside( t - t_offset_i ) * offset_i +

    post-seismic_deformation as decaying log (psd_log):   psd_eq_i * np.log( 1 + Heaviside( t - eq_i )/tau_i )

    :param model_type: string made of the key-word the parameters to be estimated.

    Key-word parameters are

    'trend','annual','semi-annual','seasonal','offset','psd_log'.

    'trend-seasonal-offset-psd_log' will do the full trajectory model.

    :param offset_dates: a list of offset_dates in decimal year

    :param eq_dates: a list of earthquake dates for which post-seismic deformation (psd_log) will be estimated

    :param H_fix: a dictionary including the name of the parameter to be hold fixed and the value.

    For instance to impose the co-seismic offset (North-East-Up) and relaxation time of 100 days for the
    first earthquake use:

    H_fix = { 'psd_log_offset_00':[10., 15., 0.] , 'psd_log_tau_00':[100., 100., 100.]}

    :param H_constraints: a dictionary including the name of the parameter to be constrained.

    For instance to impose a 50 days constraints around 500 days
    on the relaxation time of the second earthquake for all NEU components use: H_fix = { 'psd_log_tau_01':[[500.,50],
    [500.,50] , [500.,50]]}

    :param H_bounds: a dictionary including the bounds.

    For instance to impose a relaxation time for the third earthquake to be in the range
    of 2 to 3 years, for all NEU components use: H_bounds = { 'psd_log_tau_02':[[2*365.,3*365.], [[2*365.,3*365.] ,
    [[2*365.,3*365.]]}

    :param component: string , component for which the trajectory model will be estimated.

    :param verbose: verbose mode

    :note: Unlike most pyacs.gts functions, trajectory returns 4 elements: the results as a dictionary, the model Gts,
    the residual Gts and a Gts with model predictions at every day.

    """

    # import

    import numpy as np
    import pyacs.lib.astrotime as at

    # fills the H_fix, H_constraints & H_bounds for the components

    if 'N' in component:
        i = 0
        # H_fix
        H_fix_N = {}
        for k, v in H_fix.items():
            if isinstance(v, list):
                H_fix_N[k] = v[i]
            else:
                H_fix_N[k] = v
        # H_constraints
        H_constraints_N = {}
        for k, v in H_constraints.items():
            if isinstance(v[0], list):
                H_constraints_N[k] = v[i]
            else:
                H_constraints_N[k] = v
        # H_bounds
        H_bounds_N = {}
        for k, v in H_bounds.items():
            if isinstance(v[0], list):
                H_bounds_N[k] = v[i]
            else:
                H_bounds_N[k] = v

    if 'E' in component:
        i = 1
        # H_fix
        H_fix_E = {}
        for k, v in H_fix.items():
            if isinstance(v, list):
                H_fix_E[k] = v[i]
            else:
                H_fix_E[k] = v
        # H_constraints
        H_constraints_E = {}
        for k, v in H_constraints.items():
            if isinstance(v[0], list):
                H_constraints_E[k] = v[i]
            else:
                H_constraints_E[k] = v
        # H_bounds
        H_bounds_E = {}
        for k, v in H_bounds.items():
            if isinstance(v[0], list):
                H_bounds_E[k] = v[i]
            else:
                H_bounds_E[k] = v

    if 'U' in component:
        i = 2
        # H_fix
        H_fix_U = {}
        for k, v in H_fix.items():
            if isinstance(v, list):
                H_fix_U[k] = v[i]
            else:
                H_fix_U[k] = v
        # H_constraints
        H_constraints_U = {}
        for k, v in H_constraints.items():
            if isinstance(v[0], list):
                H_constraints_U[k] = v[i]
            else:
                H_constraints_U[k] = v
        # H_bounds
        H_bounds_U = {}
        for k, v in H_bounds.items():
            if isinstance(v[0], list):
                H_bounds_U[k] = v[i]
            else:
                H_bounds_U[k] = v

    # Run the estimation

    from pyacs.gts.lib.model.non_linear_gts_model import nl_gts_fit

    t_mjd = at.decyear2mjd(self.data[:, 0])
    t_mjd_ed = np.arange(t_mjd[0], t_mjd[-1])
    t_ed = at.mjd2decyear(t_mjd_ed)

    # North
    if 'N' in component:
        if verbose:
            print("-- Running trajectory model for site %s component North" % self.code)
        i = 1
        (H_res_N, model_N, residuals_N, model_ed_N) = nl_gts_fit(self.data[:, 0], self.data[:, i], self.data[:, i + 3], \
                                                                 model_type, offset_dates=offset_dates,
                                                                 eq_dates=eq_dates, \
                                                                 H_fix=H_fix_N, H_constraints=H_constraints_N,
                                                                 H_bounds=H_bounds_N, verbose=verbose)

    else:
        model_N = np.vstack((self.data[:, 0], self.data[:, 0] * 0.))
        residuals_N = np.vstack((self.data[:, 0], self.data[:, 0] * 0.))
        model_ed_N = np.vstack((t_ed, t_ed * 0.))

    # East
    if 'E' in component:
        if verbose:
            print("-- Running trajectory model for site %s component East" % self.code)
        i = 2
        (H_res_E, model_E, residuals_E, model_ed_E) = nl_gts_fit(self.data[:, 0], self.data[:, i], self.data[:, i + 3], \
                                                                 model_type, offset_dates=offset_dates,
                                                                 eq_dates=eq_dates, \
                                                                 H_fix=H_fix_E, H_constraints=H_constraints_E,
                                                                 H_bounds=H_bounds_E, verbose=verbose)

    else:
        model_E = np.vstack((self.data[:, 0], self.data[:, 0] * 0.))
        residuals_E = np.vstack((self.data[:, 0], self.data[:, 0] * 0.))
        model_ed_E = np.vstack((t_ed, t_ed * 0.))

    # Up
    if 'U' in component:
        if verbose:
            print("-- Running trajectory model for site %s component Up" % self.code)
        i = 3
        (H_res_U, model_U, residuals_U, model_ed_U) = nl_gts_fit(self.data[:, 0], self.data[:, i], self.data[:, i + 3], \
                                                                 model_type, offset_dates=offset_dates,
                                                                 eq_dates=eq_dates, \
                                                                 H_fix=H_fix_U, H_constraints=H_constraints_U,
                                                                 H_bounds=H_bounds_U, verbose=verbose)

    else:
        model_U = np.vstack((self.data[:, 0], self.data[:, 0] * 0.))
        residuals_U = np.vstack((self.data[:, 0], self.data[:, 0] * 0.))
        model_ed_U = np.vstack((t_ed, t_ed * 0.))

    # prepare return output

    H_res = [H_res_N, H_res_E, H_res_U]

    model_gts = self.copy(data_xyz=False)
    model_gts.data[:, 1] = model_N[:, -1]
    model_gts.data[:, 2] = model_E[:, -1]
    model_gts.data[:, 3] = model_U[:, -1]

    residual_gts = self.copy(data_xyz=False)
    residual_gts.data[:, 1] = residuals_N[:, -1]
    residual_gts.data[:, 2] = residuals_E[:, -1]
    residual_gts.data[:, 3] = residuals_U[:, -1]

    model_ed_gts = self.copy(data_xyz=False)
    model_ed_gts.data = np.zeros((t_ed.shape[0], 10))
    model_ed_gts.data[:, 0] = t_ed
    model_ed_gts.data[:, 1] = model_ed_N[:, -1]
    model_ed_gts.data[:, 2] = model_ed_E[:, -1]
    model_ed_gts.data[:, 3] = model_ed_U[:, -1]

    return H_res, model_gts, residual_gts, model_ed_gts
