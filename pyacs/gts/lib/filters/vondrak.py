"""
Vondrak filter for Gts.
"""


###############################################################################
def vondrak(self , fc , in_place=False , verbose=True, component='NEU'):
###############################################################################
    """
    Return a filtered Gts using a Vondrak filter.

    Parameters
    ----------
    fc : float
        Cutoff frequency in cycles per year.
    in_place : bool, optional
        If True, replace the current time series; otherwise return a new Gts.
    verbose : bool, optional
        Verbose mode.
    component : str, optional
        Components to filter ('NEU' or subset).

    Returns
    -------
    Gts
        Filtered time series.
    """


    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG
    import pyacs.debug

    ### copy
    new_gts=self.copy( data_xyz=None )

    ### filter      
    if 'N' in component:
        VERBOSE('Computing Vondrak filter for component North')
        new_gts.data[:,1] = __vondrak(self.data[:,0], self.data[:,1], fc, return_partials=False)
    
    
    if 'E' in component:
        VERBOSE('Computing Vondrak filter for component East')
        new_gts.data[:,2] = __vondrak(self.data[:,0], self.data[:,2], fc, return_partials=False)

    if 'U' in component:
        VERBOSE('Computing Vondrak filter for component Up')
        new_gts.data[:,3] = __vondrak(self.data[:,0], self.data[:,3], fc, return_partials=False)


    ### return
    if in_place:
            self = new_gts
            return self
    else:
        return  new_gts


###############################################################################
def __vondrak(t, x, fc, return_partials=False):
    """
    Pass a Vondrak filter through a 1D time series (internal routine).

    Parameters
    ----------
    t : ndarray
        Dates.
    x : ndarray
        Time series values.
    fc : float
        Cutoff frequency.
    return_partials : bool, optional
        If True, return partial derivatives dd/dy. Default is False.

    Returns
    -------
    ndarray or tuple
        Filtered time series xs; or (xs, partials) if return_partials.
    """

    import numpy
    import scipy.linalg as linalg
    import pyacs.lib.glinalg
    
    eps = (7.223147119819503*fc) ** 6 / (len(t) - 3)
    num = 6 * numpy.sqrt(t[2:-1] - t[1:-2])
    den = numpy.sqrt(t[-1] - t[0])
    
    a = numpy.hstack((0, 0, 0, num / den / ((t[0:-3] - t[1:-2]) * (t[0:-3] - t[2:-1]) * (t[0:-3] - t[3:])),   0, 0, 0))
    b = numpy.hstack((0, 0, 0, num / den / ((t[1:-2] - t[0:-3]) * (t[1:-2] - t[2:-1]) * (t[1:-2] - t[3:])),   0, 0, 0))
    c = numpy.hstack((0, 0, 0, num / den / ((t[2:-1] - t[0:-3]) * (t[2:-1] - t[1:-2]) * (t[2:-1] - t[3:])),   0, 0, 0))
    d = numpy.hstack((0, 0, 0, num / den / ((t[3:]   - t[0:-3]) * (t[3:]   - t[1:-2]) * (t[3:]   - t[2:-1])), 0, 0, 0))
    
    d0 = eps + a[3:]**2 + b[2:-1]**2 + c[1:-2]**2 + d[0:-3]**2
    d1 = a[3:-1] * b[3:-1] + b[2:-2] * c[2:-2] + c[1:-3] * d[1:-3]
    d2 = a[3:-2] * c[3:-2] + b[2:-3] * d[2:-3]
    d3 = a[3:-3] * d[3:-3]
    
    A  = numpy.diag(d0) + numpy.diag(d1,1) + numpy.diag(d1,-1) + numpy.diag(d2,2) + numpy.diag(d2,-2) + numpy.diag(d3,3) + numpy.diag(d3,-3)
    
    if (return_partials):
        A  = eps * pyacs.lib.glinalg.syminv(A)
        xs = pyacs.lib.glinalg.dot(A, x)
        return (xs, A)
    else:
        xs = linalg.solve(A, eps*x)
        return xs
