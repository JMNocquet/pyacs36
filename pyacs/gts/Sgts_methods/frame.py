###################################################################
def frame(self,frame=None,euler=None,w=None,verbose=False):
###################################################################
    """
    Rotates the time series according to an Euler pole.
    User must provide either frame, euler or w.
     
     
    :param frame: str, implemented values are 'soam','nas','nazca','inca','nas_wrt_soam','inca_wrt_soam'.
    :param euler: Euler values provided either as a \
    string 'euler_lon/euler_lat/euler_w', a list [euler_lon,euler_lat,euler_w] or \
    a 1D numpy array np.array([euler_lon,euler_lat,euler_w])
    :param w: rotation rate vector in rad/yr, provided either as a \
    string 'wx/wy/wz', a list [wx,wy,wz] or \
    a 1D numpy array np.array([wx,wy,wz])  
    
    :return: the new Sgts instance in new frame
    
    :ref: All values for frames are from Nocquet et al., Nat Geosc., 2014.
    
    """
    
    # import
    import numpy as np
    import pyacs.lib.euler
    from pyacs.gts.Sgts import Sgts
    from pyacs.gts.Gts import Gts

     
    # check arguments are OK
     
    if [frame,euler,w].count(None) != 2:
        print('!!! ERROR: define either argument frame, euler or w ')
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
        print("!!! ERROR: requested frame ",frame," not known")
        print("!!! ERROR: available frames are: ", list(lEuler.keys()))
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
            print('!!! ERROR: argument w not understood: ',w)
            return(None) 
         
        euler_vector=np.array(pyacs.lib.euler.rot2euler([w[0],w[1],w[2]]))
    
    # case euler vector
    if euler is not None:
    
        if ( isinstance(euler,str) ) and '/' in euler:
            euler=np.array(list(map(float,euler.split('/'))))
    
        if isinstance(euler,list):
            euler=np.array(euler)
         
        if not isinstance(euler,np.ndarray):
            print('!!! ERROR: argument euler not understood: ',euler)
            return(None) 
         
        euler_vector=np.array(euler)
     
    # converts the gts
     
    for gts in self.lGts():
        if verbose:print("-- Processing ",gts.code)
        try:
            new_gts=gts.remove_pole(euler_vector,pole_type='euler',in_place=False, verbose=verbose)
        except (RuntimeError, TypeError, NameError):
            print("!!! Error processing ",gts.code)
            continue
        if isinstance(new_gts,Gts):
            New_Sgts.append(new_gts)
        else:
            print("!!! Error processing ",gts.code, "!!! No time series created.")
    
    return(New_Sgts)
