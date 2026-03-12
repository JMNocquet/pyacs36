"""
Median filter for Gts based on scipy.signal.medfilt.
"""

###############################################################################
def median_filter(self , n , in_place=False , verbose=True):
###############################################################################
    """Filter time series with scipy.signal.medfilt (median filter).

    Parameters
    ----------
    n : int
        Window size (must be odd).
    in_place : bool, optional
        If True, replace self. Default is False.
    verbose : bool, optional
        Verbose mode. Default is True.

    Returns
    -------
    Gts
        Filtered time series.

    Notes
    -----
    Applied to .data; .data_xyz is set to None.
    """
    ### Check if n is odd
    if n % 2 == 0:
        raise ValueError("n must be odd")
    if n >= self.data.shape[0]:
        print("Filter length (%d) must be less than the number of data points in Gts (%d)" % (n, self.data.shape[0]))
        raise ValueError("n must be less than the number of data points")
    ### import
    import scipy.signal
    
    ### copy
    new_gts=self.copy( data_xyz=None )

    ### filter    
    new_gts.data[:,1] = scipy.signal.medfilt(self.data[:,1],n)
    new_gts.data[:,2] = scipy.signal.medfilt(self.data[:,2],n)
    new_gts.data[:,3] = scipy.signal.medfilt(self.data[:,3],n)

    ### return    
    if in_place:
            self = new_gts
            return self 
    else:
        return  new_gts 
