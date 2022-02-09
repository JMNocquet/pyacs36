"""
Wiener filter for Gts based on scipy.signal.wiener.

https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.signal.wiener.html


"""

###############################################################################
def wiener(self , in_place=False , verbose=True , my_size=15, noise=None):
###############################################################################
    """
    returns a filtered time series using scipy.signal.wiener
    
    See documentation for the filter parameters.
    
    :param in_place: if True then replace the current time series
    :param verbose: boolean, verbose mode

    :return: the filtered time series

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
