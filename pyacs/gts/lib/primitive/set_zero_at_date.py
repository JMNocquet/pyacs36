###################################################################
def set_zero_at_date(self,date, offset=None):
###################################################################
    """
    Translate the time series so that values are zero at a given date.

    If the provided date does not exist in the series, the next available date is used.

    Parameters
    ----------
    date : float
        Date in decimal year.
    offset : float or list or ndarray, optional
        Offset in mm to add. Float (same for N,E,U) or 3-element list/array for N,E,U.

    Returns
    -------
    Gts
        New Gts with translation applied; .data_xyz set to None.
    """

    # import 
    import inspect
    import numpy as np

    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG

    if self.data is None:
        ERROR("time series for site %s has no data" % self.code)
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
        ERROR("bad date provided")
        return()
    
    VERBOSE("Removing (mm) N: %8.2lf E %8.2lf U %8.2lf " % (new_data[index,1]*1E3,new_data[index,2]*1E3,new_data[index,3]*1E3))
     
    new_data[:,1]=new_data[:,1]-new_data[index,1] + cn
    new_data[:,2]=new_data[:,2]-new_data[index,2] + ce
    new_data[:,3]=new_data[:,3]-new_data[index,3] + cu

    new_Gts=self.copy(data=new_data)

    # .data_xyz set to None
    new_Gts.data_xyz = None

    return(new_Gts)
