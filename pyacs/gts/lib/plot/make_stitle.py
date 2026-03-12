###############################################################################
def make_stitle(component, info=[]):
###############################################################################
    """
    Build a subplot title string from component and optional info.

    Parameters
    ----------
    component : str
        Component to display ('N', 'E', or 'U').
    info : list or str, optional
        Additional string(s) for the subplot title.

    Returns
    -------
    str
        Title string (e.g. 'North' or custom info).
    """
    
    
    if info == '':
  
        # component
        if component=='N':
            txt_component='North'
            idx = 0
        if component=='E':
            txt_component='East'
            idx = 1
        if component=='U':
            txt_component='Up'
            idx = 2
        
        str_title = ("%s " % (txt_component) )
    
    else:

        str_title = info
    return(str_title)
