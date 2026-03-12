###################################################################
def sub(self,lexclude=[],linclude=None):
###################################################################
    """Return new Sgts with subset by exclude/include lists.

    Parameters
    ----------
    lexclude : list, optional
        Site codes to exclude. Default is [].
    linclude : list, optional
        If not None, only these codes are kept; default None means linclude is ignored.

    Returns
    -------
    Sgts
        New Sgts instance.
    """

    from pyacs.gts.Sgts import Sgts
    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG

    import inspect
    VERBOSE("Running Sgts.%s" % inspect.currentframe().f_code.co_name)

    if linclude is not None:
        sub_Sgts = Sgts(read=False)
        for code in linclude:
            if ( code not in lexclude) and ( self.has_ts( code ) ):
                sub_Sgts.append(self.__dict__[code].copy())
    else:
        sub_Sgts = self.copy()
        for code in lexclude:
            if self.has_ts( code ):
                sub_Sgts.delts(code)
     
    return( sub_Sgts )
