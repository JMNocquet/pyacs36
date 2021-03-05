###################################################################
def insert_gts_data(self,gts,in_place=False, verbose=False):
###################################################################
    """

    insert data (and/or) .data_xyz of a gts into the current gts
    
    :param gts: time series to be inserted
    :param in_place: boolean, if True add_obs to the current Gts, if False, returns a new Gts
    :param verbose: verbose mode 
  
    :return : new Gts or the modified Gts if in_place
    """
    
    # import 
    import numpy as np
    from pyacs.gts.Gts import Gts


    # check argument
    
    if not isinstance(gts,Gts) :
        raise TypeError
    
    # new_gts
    
    new_gts = self.copy()
    
    # make sur gts is ordered
    
    gts.reorder()
    
    # exclude data from the original time series
    
    new_gts = new_gts.exclude_periods([ gts.data[0,0] , gts.data[-1,0]] )

    # .data
    
    if new_gts.data is not None:
        new_gts.data = np.vstack(( new_gts.data, gts.data ))

    # .data_xyz

    if new_gts.data_xyz is not None:
        new_gts.data_xyz = np.vstack(( new_gts.data_xyz, gts.data_xyz ))

    # reorder
    
    new_gts.reorder()
    
    # in_place ?
    if in_place:
        self = new_gts
    
    return(new_gts)
