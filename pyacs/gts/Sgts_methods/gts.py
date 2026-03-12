###################################################################
def gts(self , method , *args , **kwargs):
###################################################################
    """Apply a Gts method to every series in this Sgts.

    Parameters
    ----------
    method : str
        Name of the Gts method (e.g. 'detrend', 'copy').
    *args : tuple
        Positional arguments for the method.
    **kwargs : dict
        Keyword arguments for the method.

    Returns
    -------
    Sgts
        New Sgts with results (or self for in-place methods).

    Examples
    --------
    >>> ts.gts('detrend', periods=[2010.0, 2013.0])
    """

    from pyacs.gts.Sgts import Sgts
    from tqdm import tqdm


    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG
    import pyacs.debug

    import inspect

    VERBOSE("Running Sgts.%s" % inspect.currentframe().f_code.co_name)

    verbose = kwargs.get('verbose', False)
     
    new_ts = Sgts(read=False)
     
    lsite=self.lcode()


    for site in tqdm( sorted( lsite ) , desc=method):
        try:
            func = getattr(self.__dict__[site], method)
            new_ts.append(func(*args, **kwargs) )
        except:
            WARNING("problem with method %s on gts %s. Removed from output" % ( method, site ))


    return( new_ts )
