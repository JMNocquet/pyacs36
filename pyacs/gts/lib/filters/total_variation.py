"""
Total variation filters from https://github.com/albarji/proxTV
Total variation filters are useful to preserve edges in a signal (edge filter).
"""

###############################################################################
def edge(self ,  lbda , in_place=False , verbose=True ):
###############################################################################
    """
    Edge Gts filter using a L1 total variation filter.
    The signal is assumed to be piecewise constant.
    
    :param lbda: lambda parameter
    :param in_place: if True then replace the current time series
    :param verbose: boolean, verbose mode

    :return: the filtered time series
 
    :reference: https://github.com/albarji/proxTV
 
    """
    
    import prox_tv as ptv    
    
    ### copy
    new_gts=self.copy( data_xyz=None )
   
    ### filter  
    new_gts.data[:,1] = ptv.tv1_1d(self.data[:,1]*1000.,lbda) /1000.
    new_gts.data[:,2] = ptv.tv1_1d(self.data[:,2]*1000.,lbda) /1000.
    new_gts.data[:,3] = ptv.tv1_1d(self.data[:,3]*1000.,lbda) /1000.

    ### return    
    if in_place:
            self = new_gts
            return self 
    else:
        return  new_gts 

###############################################################################
def edge_l2(self , lbda , in_place=False , verbose=True ):
###############################################################################
    """
    Gts filter using a L2 total variation filter.
    The signal is assumed to be detrended.
    
    :param lbda: lambda parameter
    :param in_place: if True then replace the current time series
    :param verbose: boolean, verbose mode

    :return: the filtered time series
    :reference: https://github.com/albarji/proxTV
 
    """
    
    import prox_tv as ptv    
    
    ### copy
    new_gts=self.copy( data_xyz=None )

    ### filter    
    
    new_gts.data[:,1] = ptv.tv2_1d(self.data[:,1]*1000.,lbda) /1000.
    new_gts.data[:,2] = ptv.tv2_1d(self.data[:,2]*1000.,lbda) /1000.
    new_gts.data[:,3] = ptv.tv2_1d(self.data[:,3]*1000.,lbda) /1000.

    ### return    
    if in_place:
            self = new_gts
            return self 
    else:
        return  new_gts 
