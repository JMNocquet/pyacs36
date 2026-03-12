###################################################################
def read_gmt(self,gmt=True,verbose=False,vel=False):
###################################################################
    """Read GMT psvelo file and set lon/lat (and optionally velocity) for each Gts.

    Parameters
    ----------
    gmt : bool or str, optional
        True = read '../stat/pyacs_void.gmt'; str = path to GMT file. Default is True.
    verbose : bool, optional
        Verbose mode. Default is False.
    vel : bool, optional
        If True, fill .velocity from file. Default is False.

    Returns
    -------
    None
        Modifies self in place.

    Notes
    -----
    Always in-place.
    """


    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG

    import inspect

    VERBOSE("Running Sgts.%s" % inspect.currentframe().f_code.co_name)

    if gmt==True:
            gmt_file=self.dir+'/../stat/pyacs_void.gmt'
    if  isinstance(gmt,str):
        gmt_file=gmt
     
    if gmt != False:
        VERBOSE("Importing modules Velocity_Field, Gpoint %s" % gmt_file)
    
        from pyacs.vel_field import Velocity_Field
         
        VERBOSE("Reading %s" % gmt_file)
        velf=Velocity_Field()
        velf.read(file_name=gmt_file)
        lsite=velf.lcode()
         
        for gts in self.lGts():
            VERBOSE("Processing %s" % gts.code)
            if gts.code in lsite:
                M=velf.site(gts.code)
                gts.lon=M.lon
                gts.lat=M.lat
                if vel:
                    gts.velocity=[M.Ve,M.Vn,0.0,M.SVe,M.SVn,0.0]
