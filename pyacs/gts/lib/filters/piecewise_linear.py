"""
Piecewise linear approximation of time series. 
Based on https://pypi.org/project/pwlf. 
"""

def pwlf(self, component , n, in_place=False, verbose=False, output=False):
    """
    Perform a piecewise approximation of a time series. Since it the routine is 1D, the component E,N, or U needs to be specified.

    :param component: component used for the decomposition. Must be 'E','N' or 'U'
    :param n: number of segments
    :param in_place: if True then replace the current time series
    :param verbose: boolean, verbose mode
    :output: if False, the predicted time series is returned. If True, then a list of dates is returned.

    """

    # import
    import pwlf
    
    ### copy
    new_gts=self.copy( data_xyz=None )
    
    # select component
    
    ic = 0
    if component == 'E':
        ic = 2
    if component == 'N':
        ic = 1
    if component == 'U':
        ic = 3
    
    # initialize
    my_pwlf = pwlf.PiecewiseLinFit(self.data[:,0], self.data[:,ic]) 
    
    #perform fit and get the dates
    np_dates = my_pwlf.fit(n)
    
    # return of output is True
    if output: 
        return np_dates
    
    # otherwise just continue
    # make prediction for the fitted component
    
    new_gts.data[:,ic] = my_pwlf.predict( new_gts.data[:,0] )
    
    # make the prediction for the other components
    
    # uncertainties to 1 mm
    
    new_gts.data[:,4:7] = 1. * 1.E-3
    new_gts.data[:,7:]  = 0.
    
    for ic in [1,2,3]:
    
        my_pwlf = pwlf.PiecewiseLinFit(self.data[:,0], self.data[:,ic]) 
        my_pwlf.fit_with_breaks( np_dates )    
        new_gts.data[:,ic] = my_pwlf.predict( new_gts.data[:,0] )
    

    ### return    
    if in_place:
            self = new_gts
            return self 
    else:
        return  new_gts 
