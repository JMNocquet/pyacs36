"""
Minimum component filter for Gts.

"""
###############################################################################
def minimum_component(self , mask_period=[], p=1, fcut=None, Q=None , in_place=False , verbose=True ):
###############################################################################
    """
    Minimum component filtering for Gts.
    Minimum component filtering is useful for determining the background
    component of a signal in the presence of spikes
    :param mask_periods: periods (list or list of lists) which should be ignored for smoothing
    :param p: integer (optional). polynomial degree to be used for the fit (default = 1)
    :param fcut: float (optional). the cutoff frequency for the low-pass filter.  Default value is f_nyq / sqrt(N)
    :param Q: float (optional). the strength of the low-pass filter.  Larger Q means a steeper cutoff. default value is 0.1 * fcut
    :param in_place: if True then replace the current time series
    :param verbose: boolean, verbose mode
    :return: the filtered time series
    :note:
    This code follows the procedure explained in the book
    "Practical Statistics for Astronomers" by Wall & Jenkins book, as
    well as in Wall, J, A&A 122:371, 1997
    """
    
    # import
    import numpy as np
    
    # handle masked periods

    ###############################################################################
    def __ensure_list_of_list(ll):
    ###############################################################################
        """
        Ensures ll is a list of lists
        
        [a,b] returns [[a,b]], and [[a,b]] returns [[a,b]]
        """
    
        # check this is a list
        
        if not isinstance(ll,list):
            raise TypeError('!!! __ensure_list_of_list requires a list or a list of list as argument: ',ll)
        
        if ll == []:
            return [[]]
        
        # case simple list
        if not isinstance(ll[0],list):
            return([ll])
        # case list of list
        else:
            return(ll)
    ###############################################################################
    
    # ensure lperiod is a list of lists
    lperiod = __ensure_list_of_list( mask_period )
  
    # ma = feature mask 
    feature_mask = np.copy( self.data[:,0] )
    feature_mask[:] = False

    np_index=np.array([],dtype=int)
    
    # check lperiod is not [[]]
    if lperiod != [[]]:
        # masked periods
        for period in lperiod:
            
            # make actual extraction - case data_xyz now handled
            start_date_period= period[0]
            end_date_period  = period[1]
    
            np_idx = np.where((self.data[:,0]>=start_date_period) & (self.data[:,0] <=end_date_period ))
            
            np_index = np.append( np_index , np_idx )

        feature_mask[ np_index ] = True

    ### copy
    new_gts=self.copy( data_xyz=None )


    ### run filter 
    new_gts.data[:,1] = __min_component_filter(self.data[:,0],self.data[:,1], feature_mask , \
                                                   p=p, \
                                                   fcut=fcut , \
                                                   Q=Q)
    
    
    new_gts.data[:,2] = __min_component_filter(self.data[:,0],self.data[:,2], feature_mask , \
                                                   p=p, \
                                                   fcut=fcut , \
                                                   Q=Q)

    new_gts.data[:,3] = __min_component_filter(self.data[:,0],self.data[:,3], feature_mask , \
                                                   p=p, \
                                                   fcut=fcut , \
                                                   Q=Q)


    # restore original values for masked periods

    new_gts.data[np_index,1] = np.copy( self.data[np_index,1] )
    new_gts.data[np_index,2] = np.copy( self.data[np_index,2] )
    new_gts.data[np_index,3] = np.copy( self.data[np_index,3] )
    

    ### return
    if in_place:
            self = new_gts
            return self
    else:
        return  new_gts

###############################################################################
def __min_component_filter(x, y, feature_mask, p=1, fcut=None, Q=None):
###############################################################################

    """
    Minimum component filtering

    Minimum component filtering is useful for determining the background
    component of a signal in the presence of spikes

    Parameters
    x : array_like
        1D array of evenly spaced x values
    y : array_like
        1D array of y values corresponding to x
    feature_mask : array_like, same size as x & y 
        1D mask array giving the locations of features in the data which
        should be ignored for smoothing
    p : integer (optional)
        polynomial degree to be used for the fit (default = 1)
    fcut : float (optional)
        the cutoff frequency for the low-pass filter.  Default value is
        f_nyq / sqrt(N)
    Q : float (optional)
        the strength of the low-pass filter.  Larger Q means a steeper cutoff
        default value is 0.1 * fcut

    Returns
    y_filtered : ndarray
        The filtered version of y.

    Notes
    This code follows the procedure explained in the book
    "Practical Statistics for Astronomers" by Wall & Jenkins book, as
    well as in Wall, J, A&A 122:371, 1997
    """
    
    import numpy as np
    from scipy import fftpack
    
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    
    feature_mask = np.asarray(feature_mask, dtype=bool)

    if ((x.ndim != 1) or (x.shape != y.shape) or (y.shape !=
                                                  feature_mask.shape)):
        raise ValueError('x, y, and feature_mask must be 1 dimensional '
                         'with matching lengths')

    if fcut is None:
        f_nyquist = 1. / (x[1] - x[0])
        fcut = f_nyquist / np.sqrt(len(x))

    if Q is None:
        Q = 0.1 * fcut

    # compute polynomial features
    XX = x[:, None] ** np.arange(p + 1)

    # compute least-squares fit to non-masked data
    beta = np.linalg.lstsq(XX[~feature_mask], y[~feature_mask])[0]

    # subtract polynomial fit and mask the data
    y_mask = y - np.dot(XX, beta)
    y_mask[feature_mask] = 0

    # get Fourier transforms of arrays
    yFT_mask = fftpack.fft(y_mask)

    # compute (shifted) frequency array for filter
    N = len(x)
    f = fftpack.ifftshift((np.arange(N) - N / 2.) * 1. / N / (x[1] - x[0]))

    # construct low-pass filter
    filt = np.exp(- (Q * (abs(f) - fcut) / fcut) ** 2)
    filt[abs(f) < fcut] = 1

    # reconstruct filtered signal
    y_filtered = fftpack.ifft(yFT_mask * filt).real + np.dot(XX, beta)

    return y_filtered


