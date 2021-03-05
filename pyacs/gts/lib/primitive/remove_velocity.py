###################################################################
def remove_velocity(self,vel_neu,in_place=False):
###################################################################
    """
    remove velocity from a time series
    vel_neu is a 1D array of any arbitrary length, but with the velocities (NEU) to be removed in the first 3 columns
    if in_place = True then replace the current time series
    """

    # import 
    import inspect
    import numpy as np

    # check data is not None
    from pyacs.gts.lib.errors import GtsInputDataNone
    
    try:
        if self.data is None:
            # raise exception
            raise GtsInputDataNone(inspect.stack()[0][3],__name__,self)
    except GtsInputDataNone as error:
        # print PYACS WARNING
        print( error )
        return( self )

    if self.t0 != None:
        ref_date=self.t0
    else:    
        ref_date=np.mean(self.data[0,0])
    
    v_neu= np.array( vel_neu[0:3] )
    
    
    new_data=self.data.copy()
     
    new_data[:,1:4]=new_data[:,1:4]-np.dot( (new_data[:,0]-ref_date).reshape(-1,1) , v_neu.reshape(1,3) )

    new_Gts=self.copy(data=new_data)
    
    # .data_xyz set to None
    new_Gts.data_xyz = None

    vel = vel_neu
    
    new_Gts.velocity = vel
    
    if in_place:
        self.data=new_Gts.data.copy()
        return( self )
    else:
        return(new_Gts)
