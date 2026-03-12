def ivel(self):
    """Compute instantaneous velocity (time series should be filtered first).

    Returns
    -------
    Gts
        New Gts instance with velocity; units m/yr.

    Notes
    -----
    Output dates are shifted by 1 day so the first ivel date is the
    second date of the input time series, and so on.
    """
    
    import numpy as np
    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG
    import pyacs.debug

    fts = self

    # compute the derivative
    x = fts.data[:,0]
    dx = x[:-1] + np.diff(x)
    # north
    fy = fts.data[:,1]*1.E3
    dyn = np.diff(fy)/(np.diff(x))
    # east
    fy = fts.data[:,2]*1.E3
    dye = np.diff(fy)/(np.diff(x))
    # up
    fy = fts.data[:,3]*1.E3
    dyu = np.diff(fy)/(np.diff(x))

    new_gts = fts.copy()
    new_gts.data = np.zeros((dx.shape[0],10))+1.E-3
    new_gts.data[:,0] = dx
    new_gts.data[:,1] = dyn*1E-3
    new_gts.data[:,2] = dye*1E-3
    new_gts.data[:,3] = dyu*1E-3

    return new_gts