###################################################################
def gts(self , method , *args , **kwargs):
###################################################################
    """
    apply a gts method to all Gts instance of the current Sgts object
     
    :param method: Gts method to be applied as string
    :param *arg: arguments for the Gts method to be applied 
    :param **kwarg: keyword arguments for the Gts method to be applied 
    
    :example : ts.gts('detrend',periods=[2010.0,2013.0])
    """

    from pyacs.gts.Sgts import Sgts

    verbose = kwargs.get('verbose', False)        
     
    new_ts = Sgts(read=False)
     
    lsite=self.lcode()
     
    for site in sorted( lsite ):
        if verbose:
            print('-- processing ', site)
        
        try:
            func = getattr(self.__dict__[site], method)
            new_ts.append(func(*args, **kwargs) )
        except:
            print('!!!WARNING: problem with method %s on gts %s. Removed from output' % ( method, site ))
    return( new_ts )
