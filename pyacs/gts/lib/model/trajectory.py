"""
Non-linear trajectory models for geodetic time series.
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
    Estimate parameters of a (non-linear) trajectory model for the time series.

    Model: y(t) = trend + annual + semi-annual + offsets + post-seismic (psd_log).

    Parameters
    ----------
    model_type : str
        String of key-word parameters to estimate: 'trend', 'annual', 'semi-annual',
        'seasonal', 'offset', 'psd_log'. E.g. 'trend-seasonal-offset-psd_log' for full model.
    offset_dates : list, optional
        List of offset dates in decimal year.
    eq_dates : list, optional
        List of earthquake dates for post-seismic (psd_log) estimation.
    H_fix : dict, optional
        Parameters to fix, e.g. {'psd_log_offset_00': [10., 15., 0.], 'psd_log_tau_00': [100., 100., 100.]}.
    H_constraints : dict, optional
        Parameters to constrain (center, sigma), e.g. {'psd_log_tau_01': [[500., 50], [500., 50], [500., 50]]}.
    H_bounds : dict, optional
        Bounds, e.g. {'psd_log_tau_02': [[2*365., 3*365.], ...]}.
    component : str, optional
        Components to estimate ('NEU' or subset).
    verbose : bool, optional
        Verbose mode.

    Returns
    -------
    tuple
        (results_dict, model_Gts, residual_Gts, daily_predictions_Gts). Unlike most pyacs.gts
        functions, trajectory returns these 4 elements.
    """

    # import

    import numpy as np
    import pyacs.lib.astrotime as at
    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG


    # after this method .data  and .data_xyz are not consistent so .data_xyz is set to None
    #self.data_xyz = None

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
        VERBOSE("Running trajectory model for site %s component North: %s" % (self.code, 'N'))
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
        VERBOSE("Running trajectory model for site %s component East: %s" % (self.code, 'E'))
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
        VERBOSE("Running trajectory model for site %s component Up: %s" % (self.code, 'U'))
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
