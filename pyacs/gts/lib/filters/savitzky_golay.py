"""
Savitzky-Golay filter for Gts based on scipy.signal.medfilt.
http://https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.savgol_filter.html#scipy.signal.savgol_filter
Additional information on Savitzky-Golay filter from https://scipy-cookbook.readthedocs.io/items/SavitzkyGolay.html

The Savitzky Golay filter is a particular type of low-pass filter, well adapted for data smoothing. For further information see: http://www.wire.tu-bs.de/OLDWEB/mameyer/cmr/savgol.pdf (or http://www.dalkescientific.com/writings/NBN/data/savitzky_golay.py for a pre-numpy implementation).

It has the advantage of preserving the original shape and
features of the signal better than other types of filtering
approaches, such as moving averages techniques.
    
The Savitzky-Golay is a type of low-pass filter, particularly
suited for smoothing noisy data. The main idea behind this
approach is to make for each point a least-square fit with a
polynomial of high order over a odd-sized window centered at
the point.
"""

###############################################################################
def savitzky_golay(self , in_place=False , verbose=True , window_length=15, polyorder=3, deriv=0, delta=1.0, mode='interp', cval=0.0):
###############################################################################
    """
    returns a filtered time series using scipy.signal.savgol_filter
    
    See documentation for the filter parameters.
    http://https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.savgol_filter.html#scipy.signal.savgol_filter
    
    :param in_place: if True then replace the current time series
    :param verbose: boolean, verbose mode

    :return: the filtered time series

    """
    
    import scipy.signal
    
    ### copy
    new_gts=self.copy( data_xyz=None )
    
    new_gts.data[:,1] = scipy.signal.savgol_filter(self.data[:,1], \
                                                   window_length=window_length, \
                                                   polyorder=polyorder, \
                                                   deriv=deriv, \
                                                   delta=delta, \
                                                   axis=-1, \
                                                   mode=mode, \
                                                   cval=cval)
    
    
    new_gts.data[:,2] = scipy.signal.savgol_filter(self.data[:,2],\
                                                    window_length=window_length, \
                                                   polyorder=polyorder, \
                                                   deriv=deriv, \
                                                   delta=delta, \
                                                   axis=-1, \
                                                   mode=mode, \
                                                   cval=cval)

    new_gts.data[:,3] = scipy.signal.savgol_filter(self.data[:,3],\
                                                   window_length=window_length, \
                                                   polyorder=polyorder, \
                                                   deriv=deriv, \
                                                   delta=delta, \
                                                   axis=-1, \
                                                   mode=mode, \
                                                   cval=cval)


    ### return    
    if in_place:
            self = new_gts
            return self 
    else:
        return  new_gts 
