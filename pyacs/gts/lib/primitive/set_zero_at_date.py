###################################################################
def set_zero_at_date(self,date,offset=None,in_place=False):
###################################################################
    """
    make a translation of a time series, setting to 0 at a given date
    if the provided date does not exist, uses the next date available
    
    :param date: date in decimal year
    :param offset: an offset (in mm) to be added. Could be a float, a list or 1D numpy array with 3 elements 
    
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

    if offset is None:
        [cn,ce,cu] =[0.0,0.0,0.0]
    else:
        if isinstance(offset,float):
            [cn,ce,cu] =[offset/1000.,offset/1000.,offset/1000.]
        elif isinstance(offset,list):
            [cn,ce,cu] =[offset[0]/1000.,offset[1]/1000.,offset[2]/1000.]
        elif isinstance(offset,np.ndarray):
            [cn,ce,cu] =[offset[0]/1000.,offset[1]/1000.,offset[2]/1000.]
    
    
    
    
    new_data=self.data.copy()

    lindex=np.where(new_data[:,0] >= date)

    try:
        index=lindex[0][0]
    except:
        print("!!! bad date provided")
        return()
    
    print("-- Removing ",new_data[index,1],new_data[index,2],new_data[index,3])
     
    new_data[:,1]=new_data[:,1]-new_data[index,1] + cn
    new_data[:,2]=new_data[:,2]-new_data[index,2] + ce
    new_data[:,3]=new_data[:,3]-new_data[index,3] + cu

    new_Gts=self.copy(data=new_data)

    # .data_xyz set to None
    new_Gts.data_xyz = None
    
    # origin XYZ would need to be changed
    
    if in_place:
        self.data=new_Gts.data.copy()
    return(new_Gts)
