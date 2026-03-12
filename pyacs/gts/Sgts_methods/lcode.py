###################################################################
def lcode(self,lexclude=[],linclude=[],min_obs=None,max_obs=None, min_duration=None, force_4_digit = False):
###################################################################
    """Return list of Gts codes in this Sgts (with optional filters).

    Parameters
    ----------
    lexclude : list, optional
        Site codes to exclude. Default is [].
    linclude : list, optional
        If non-empty, only these codes are included. Default is [].
    min_obs : int, optional
        Minimum number of epochs per site. Default is None.
    max_obs : int, optional
        Maximum number of epochs per site. Default is None.
    min_duration : float, optional
        Minimum duration in decimal year. Default is None.
    force_4_digit : bool, optional
        If True, keep only 4-character code. Default is False.

    Returns
    -------
    list
        Site codes.
    """
    # import
    import re

    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG

    lcode=[]
    if linclude != []:
        for k in list(self.__dict__.keys()):
            if ((re.match('[A-Z0-9_]{4}',k))) and (k not in lexclude) and (k in linclude):lcode.append(k)
    else:
        for k in list(self.__dict__.keys()):
            if ((re.match('[A-Z0-9_]{4}',k))) and (k not in lexclude):lcode.append(k)

    # min_obs
    if min_obs is not None:
        llcode = []
        for code in lcode:
            if self.__dict__[code].data.shape[0] < min_obs:
                VERBOSE("%s has less than %d epochs. Excluding it" %(code,self.__dict__[code].data.shape[0]))
            else:
                llcode.append(code)
        lcode = llcode

    # max_obs
    if max_obs is not None:
        llcode = []
        for code in lcode:
            if self.__dict__[code].data.shape[0] <= max_obs:
                llcode.append(code)
        lcode = llcode

    # min_duration
    if min_duration is not None:
        llcode = []
        for code in lcode:
            duration = self.__dict__[code].data.shape[-1,0] - self.__dict__[code].data.shape[0,0]
            if duration < min_duration:
                VERBOSE("%s time series has a duration less than %.3lf years. Excluding it" % (code, duration))
            else:
                llcode.append(code)
        lcode = llcode

    # 4digit case

    if force_4_digit:
        lcode = [x[:4] for x in lcode]

    return(sorted(lcode))
