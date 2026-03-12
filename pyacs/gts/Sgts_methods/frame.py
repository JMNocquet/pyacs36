###################################################################
def frame(self,frame=None,euler=None,w=None,verbose=False):
###################################################################
    """Rotate time series according to an Euler pole.

    Provide exactly one of: frame, euler, or w.

    Parameters
    ----------
    frame : str, optional
        Named frame: 'soam', 'nas', 'nazca', 'inca', 'nas_wrt_soam', 'inca_wrt_soam'.
    euler : str or list or numpy.ndarray, optional
        Euler pole: string 'lon/lat/w', list [lon,lat,w], or 1D array.
    w : str or list or numpy.ndarray, optional
        Rotation rate (rad/yr): string 'wx/wy/wz', list [wx,wy,wz], or 1D array.

    Returns
    -------
    Sgts
        New Sgts in the chosen frame.

    Notes
    -----
    Frame values from Nocquet et al., Nat. Geosci., 2014.
    """


    # import
    import numpy as np
    import pyacs.lib.euler
    from pyacs.gts.Sgts import Sgts
    from pyacs.gts.Gts import Gts

    from icecream import ic
    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG


    import inspect
    VERBOSE("Running Sgts.%s" % inspect.currentframe().f_code.co_name)

    # check arguments are OK
     
    if [frame,euler,w].count(None) != 2:
        ERROR("argument frame must be euler or w ")
        return(None)
     
    # Euler poles taken from pygvel_pole_info.py
    lEuler={}
    # from Jarrin et al., 2022
    lEuler['soam']=[ -133.29 , -18.41 , 0.121]
    lEuler['nas']=[-133.999, 12.404, 0.158]
    # from Nocquet et al., 2021
    lEuler['soam_nocquet_2014']=[-132.21,-18.83,0.121]
    lEuler['nas_nocquet_2014']=[-97.52,6.48,0.359]
    lEuler['nazca_nocquet_2014']=[-94.4,61.0,0.57]
    lEuler['inca_nocquet_2014']=[-103.729,-1.344,0.1659]
    lEuler['nas_wrt_soam_nocquet_2014']=[-83.40,15.21,0.287]
    lEuler['inca_wrt_soam_nocquet_2014']=[-63.76,22.47,0.092]
     
    # check frame case is OK
    if ( frame not in list(lEuler.keys())) and ( frame is not None):
        ERROR("requested frame %s not known" % frame)
        WARNING("available frames are: %s" % " ".join(list(lEuler.keys())))
        return(None)
     
    # initialize new gts
     
    New_Sgts=Sgts(read=False)
    
    # convert to Euler vector whatever the provided argument
     
    # case frame
    if frame is not None:
        euler_vector=np.array(lEuler[frame])
    
    # case w as rotation rate vector
    if w != None:
    
        if ( isinstance(w,str) ) and '/' in w:
            w=np.array(list(map(float,w.split('/'))))
    
        if isinstance(w,list):
            w=np.array(w)
         
        if not isinstance(w,np.ndarray):
            ERROR("argument w not understood: %s" % str(w))
            return(None) 
         
        euler_vector=np.array(pyacs.lib.euler.rot2euler(w[0],w[1],w[2]))
    
    # case euler vector
    if euler is not None:
    
        if ( isinstance(euler,str) ) and '/' in euler:
            euler=np.array(list(map(float,euler.split('/'))))
    
        if isinstance(euler,list):
            euler=np.array(euler)
         
        if not isinstance(euler,np.ndarray):
            ERROR("argument euler not understood: %s" % str(euler))
            return(None) 
         
        euler_vector=np.array(euler)
     
    # converts the gts
     
    for gts in self.lGts():
        VERBOSE("Processing %s"% gts.code)

        try:
            new_gts=gts.remove_pole(euler_vector,pole_type='euler', verbose=verbose)
        except (RuntimeError, TypeError, NameError):
            ERROR("processing %s "% gts.code)
            continue
        if isinstance(new_gts,Gts):
            New_Sgts.append(new_gts)
        else:
            ERROR("processing %s. No time series created." % gts.code)

    return(New_Sgts)
