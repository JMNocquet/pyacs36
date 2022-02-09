###############################################################################
def make_stitle(component, info=[]):
###############################################################################
    """
    
    :param component: component that will be displayed as subplottitle
    :param info: string be displayed as subplottitle
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
