def to_pytrf(self):
    """
    Convert pyacs Gts to pytrf ts object.

    Uses the pytrf library from https://github.com/prebischung/pytrf.


    Returns
    -------
    pytrf.ts.ts
        Time series in pytrf format (East, North, Up; MJD; covariances).
    """
    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG


    try:
        from pytrf.ts import ts
    except ImportError:
        ERROR("pytrf is not installed, please install it to use this function", exit=True)
    import pyacs.lib.astrotime as at
    import numpy as np
    
    # Create empty pytrf ts object
    r = ts()

    # date
    t = at.decyear2mjd(self.data[:, 0])

    # position - pytrf ts expects (east, north, up)
    y = self.data[:, [2,1,3]]
    (n, nd) = y.shape
    # covariances
    Q = np.zeros((n, nd, nd))
    Q[:, 0, 0] = self.data[:, 5]**2
    Q[:, 1, 1] = self.data[:, 4]**2
    Q[:, 2, 2] = self.data[:, 6]**2
    # covariances
    Q[:, 0,1 ] = self.data[:, 7] * np.sqrt(Q[:, 0, 0] * Q[:, 1, 1])
    Q[:, 0,2 ] = self.data[:, 9] * np.sqrt(Q[:, 0, 0] * Q[:, 2, 2])
    Q[:, 1,2 ] = self.data[:, 8] * np.sqrt(Q[:, 1, 1] * Q[:, 2, 2])
    # symmetrize
    Q[:, 1, 0] = Q[:, 0, 1]
    Q[:, 2, 0] = Q[:, 0, 2]
    Q[:, 2, 1] = Q[:, 1, 2]

    # Create ts instance
    r = ts(y, t=t, T=None, Q=Q, tunit='d', yunit='m', dims=('East', 'North', 'Up'), rotate=False, t0=None, dtrd=None)
            
    return r
