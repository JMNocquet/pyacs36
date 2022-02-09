"""
Median filter for Gts based on scipy.signal.medfilt.
"""

###############################################################################
def median_filter(self , n , in_place=False , verbose=True):
###############################################################################
    """
    returns a filtered time series using scipy.signal.medfilt
    
    :param n: size of the median filter window (must be odd)
    :param in_place: if True then replace the current time series
    :param verbose: boolean, verbose mode
    :return: the filtered time series

    :note: the filter is applied to .data and .data_xyz is set to None
    """
    
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
