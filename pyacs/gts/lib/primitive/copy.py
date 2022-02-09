###################################################################
def copy(self,data=True, data_xyz=True, loutliers=True):
###################################################################
    """
    makes a (deep) copy of the time series.
    
    By default, all attributes are also copied, including .data, .data_xyz, loutliers etc.

    Default behaviour can be modified for the following attribute:
    
    :param data: can be set to None or a 2D numpy array of shape (n,10)
    :param data_xyz: can be set to None or a 2D numpy array of shape (n,10)
    :param loutliers: False will not copy the loutliers atrribute
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

    if isinstance(new_Gts.data,np.ndarray):

        ldate_outliers = self.data[:,0][self.outliers]
        lupdated_outliers = pyacs.gts.Gts.get_index_from_dates(ldate_outliers, new_Gts.data, tol=0.01)

        new_Gts.outliers = lupdated_outliers 

    else:
        new_Gts.outliers = []

    return( new_Gts )
