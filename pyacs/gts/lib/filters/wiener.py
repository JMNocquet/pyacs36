"""
Wiener filter for Gts based on scipy.signal.wiener.

https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.signal.wiener.html


"""

###############################################################################
def wiener(self , in_place=False , verbose=True , my_size=15, noise=None):
###############################################################################
    """
    Return a filtered time series using scipy.signal.wiener.

    Parameters
    ----------
    in_place : bool, optional
        If True, replace the current time series; otherwise return a new Gts.
    verbose : bool, optional
        Verbose mode.
    my_size : int, optional
        Size of the Wiener filter window.
    noise : float, optional
        Noise estimate for the filter.

    Returns
    -------
    Gts
        Filtered time series.

    See Also
    --------
    scipy.signal.wiener
    """
    
    import scipy.signal
    
    new_gts=self.copy()
    
    new_gts.data[:,1] = scipy.signal.wiener(self.data[:,1], \
                                                   mysize=my_size, \
                                                   noise=noise)
    
    
    new_gts.data[:,2] = scipy.signal.wiener(self.data[:,2],\
                                                   mysize=my_size, \
                                                   noise=noise)

    new_gts.data[:,3] = scipy.signal.wiener(self.data[:,3],\
                                                   mysize=my_size, \
                                                   noise=noise)


    if in_place:
            self.data=new_gts.data
    else:
        return(new_gts)
