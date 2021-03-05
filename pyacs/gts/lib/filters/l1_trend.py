"""
L1 trend filter
    The L1 trend filtering method produces trend estimates that are piecewise linear from the time series 𝑦.
    :credit: https://www.cvxpy.org/examples/applications/l1_trend_filter.html
    :reference: http://stanford.edu/~boyd/papers/l1_trend_filter.html
"""


###############################################################################
def l1_trend(self , vlambda , gap=10, in_place=False , verbose=True, component='NEU'):
###############################################################################
    """
    return a piecewise linear filtered Gts
    
    :param vlambda: weight of regularization
    :param gap: gap in days to split the time series before applying the filter.default is 10.
    :param in_place: if True then replace the current time series
    :param verbose: boolean, verbose mode
    :param component: string. Default 'NEU'
    :return: the filtered time series
    :note: if there are less than 4 points in a segment, return an L1 estimated trend
    
    """
    # TODO
    # case a segment has obly 1 data is not handled yet

    # import
    from pyacs.gts.lib.message import message
    import numpy as np
    import pyacs.lib.astrotime as at


    # split time series according to gap
    lgts = self.split_gap(gap=gap, verbose=verbose)

    ### copy
    new_gts=self.copy( data_xyz=None )
    new_gts.data_xyz = None
    new_gts.data = None


    if verbose:
        print("--- time series is be filtered for %d segments" % len(lgts))

    for wgts in lgts:

        if verbose:
            datetime1 = at.decyear2datetime( wgts.data[0,0] )
            datetime2 = at.decyear2datetime( wgts.data[-1,0] )

            print("---- period %s -- %s " % (datetime1.isoformat(), \
                                             datetime2.isoformat()))



        new_data = np.zeros((wgts.data.shape[0],10))
        new_data[:,0] = np.copy( wgts.data[:,0] )
        new_data[:,4:7] = 1.E-3

        try:
            ### filter
            if 'N' in component:
                if verbose:
                    message('Computing L1 trend filter for component North')
                new_data[:,1] = __l1_trend(wgts.data[:,0], wgts.data[:,1], vlambda, verbose=verbose)
            else:
                new_data[:, 1] = wgts.data[:,1]

            if 'E' in component:
                if verbose:
                    message('Computing L1 trend filter for component East')
                new_data[:, 2] = __l1_trend(wgts.data[:, 0], wgts.data[:, 2], vlambda, verbose=verbose)
            else:
                new_data[:, 2] = wgts.data[:,2]


            if 'U' in component:
                if verbose:
                    message('Computing L1 trend filter for component Up')
                new_data[:, 3] = __l1_trend(wgts.data[:, 0], wgts.data[:, 3], vlambda, verbose=verbose)
            else:
                new_data[:, 3] = wgts.data[:,3]

        except:
            # L1 trend if less than 4 data
            if verbose:
                message('not enough data for this segment. Doing an L1 detrend')
            l1 = wgts.detrend(method='L1')
            if 'N' in component:
                bias = np.mean( l1.data[:, 1])
                new_data[:, 1] = wgts.data[:, 1] - l1.data[:, 1] + bias
            if 'E' in component:
                bias = np.mean( l1.data[:, 2])
                new_data[:, 2] = wgts.data[:, 2] - l1.data[:, 2] + bias
            if 'U' in component:
                bias = np.mean( l1.data[:, 3])
                new_data[:, 3] = wgts.data[:, 3] - l1.data[:, 3] + bias


        # adds interpolation to the time series
        if new_gts.data is None:
            new_gts.data = new_data
        else:
            new_gts.data = np.vstack((new_gts.data,new_data))


### return
    if in_place:
            self = new_gts
            return self
    else:
        return  new_gts


###############################################################################
def __l1_trend(t, y, vlambda, verbose=False):
    """

    :param t: dates an 1-D numpy array
    :param y: time series value as 1-D numpy array
    :param vlambda: weight for filtering
    :param verbose: boulean verbose mode
    :return: filtered time series
    """

    # import
    import l1tf
    import numpy as np

    return l1tf.l1_filter(np.ascontiguousarray(y),vlambda)


###############################################################################
def __l1_trend_cvxopt(t, y, vlambda, verbose=False):
    """

    :param t: dates an 1-D numpy array
    :param y: time series value as 1-D numpy array
    :param vlambda: weight for filtering
    :param verbose: boulean verbose mode
    :return: filtered time series
    """

    # import
    import numpy as np
    import cvxpy as cp
    import scipy as scipy
    import cvxopt as cvxopt
    from pyacs.gts.lib.message import message

    n = y.size

    # Form second difference matrix.
    e = np.ones((1, n))
    D = scipy.sparse.spdiags(np.vstack((e, -2*e, e)), range(3), n-2, n)

    # Solve l1 trend filtering problem.
    x = cp.Variable(shape=n)
    obj = cp.Minimize(0.5 * cp.sum_squares(y - x)
                      + vlambda * cp.norm(D*x, 1) )
    prob = cp.Problem(obj)

    # ECOS and SCS solvers fail to converge before
    # the iteration limit. Use CVXOPT instead.
    if verbose:
        message( 'solving L1 trend')

    prob.solve(solver=cp.CVXOPT, verbose=True)

    if verbose:
        my_message = 'Solver status: {}'.format(prob.status)
        message( my_message )

    # Check for error.
    if prob.status != cp.OPTIMAL:
        raise Exception("Solver did not converge!")

    if verbose:
        my_message = "optimal objective value: {}".format(obj.value)
        message( my_message )

    return np.array( x.value )