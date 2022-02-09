###################################################################
def rotate(self,angle,in_place=False):
###################################################################
    """
    rotates the axis by an angle

    :param angle: angle in decimal degrees clockwise
    
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

    # computes the rotates components
    angle_radian = np.radians( -angle )
    
    new_data=self.data.copy()
    
    new_x = new_data[:,2] * np.cos( angle_radian ) - new_data[:,1] * np.sin( angle_radian ) 
    new_y = new_data[:,2] * np.sin( angle_radian ) + new_data[:,1] * np.cos( angle_radian ) 
     
    new_data[:,1] = new_y
    new_data[:,2] = new_x

    new_data[:,4] = new_data[:,5] = 1.E-3
    new_data[:,7:] = 0.0

    new_Gts=self.copy(data=new_data)

    # .data_xyz set to None
    new_Gts.data_xyz = None
    
    
    if in_place:
        self.data=new_Gts.data.copy()
        return( self )
    else:
        return(new_Gts)
