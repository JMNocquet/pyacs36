###################################################################
def copy(self,data=True, data_xyz=True, loutliers=True):
###################################################################
    """Return a deep copy of the time series.

    By default copies .data, .data_xyz, outliers, etc. Behaviour can be
    overridden per attribute below.

    Parameters
    ----------
    data : bool or numpy.ndarray, optional
        True = copy .data; False = set to None; or (n,10) array. Default is True.
    data_xyz : bool or numpy.ndarray, optional
        True = copy .data_xyz; False = None; or (n,10) array. Default is True.
    loutliers : bool, optional
        If False, do not copy outliers. Default is True.

    Returns
    -------
    Gts
        New Gts instance.
    """

    ### import 
    
    import numpy as np
    import copy
    import pyacs.gts
    
    ### deep copy
    new_Gts=copy.deepcopy(self)
    
    ### handles .data
    
    if isinstance( data , bool):
        if not data:
            new_Gts.data = None
    else:
        new_Gts.data = data
    
    ### handles .data_xyz
    
    if isinstance( data_xyz , bool):
        if not data_xyz:
            new_Gts.data_xyz = None
    else:
        new_Gts.data_xyz = data_xyz
    
    ### handles outliers

    if isinstance(new_Gts.data,np.ndarray) and loutliers:

        ldate_outliers = self.data[:,0][self.outliers]
        lupdated_outliers = pyacs.gts.Gts.get_index_from_dates(ldate_outliers, new_Gts.data, tol=0.01)

        new_Gts.outliers = lupdated_outliers 

    else:
        new_Gts.outliers = []

    ### handles other attributes

    lattributes = ['offsets_dates', 'offsets_values', 'annual', 'semi_annual', 'velocity']
    for attribute in lattributes:
        if hasattr(self, attribute):
            setattr(new_Gts, attribute, getattr(self, attribute))

    return( new_Gts )
