###################################################################
def medvel(self,outdir=None,verbose=False):
###################################################################
    """Automatic velocity estimates using median estimator (MIDAS-style, Blewitt et al., 2016).

    Fills the velocity attribute of every Gts. Optionally writes results to outdir.

    Parameters
    ----------
    outdir : str, optional
        Output directory for results. Default is None.
    verbose : bool, optional
        Verbose mode. Default is False.

    Returns
    -------
    Sgts
        Modified Sgts (velocity set for each Gts).

    References
    ----------
    Blewitt et al. (2016), J. Geophys. Res. Solid Earth, 121(3), 2054-2068.
    """

    # import    
    from pyacs.gts.Sgts import Sgts
    import sys , os

    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG



    import inspect

    VERBOSE("Running Sgts.%s" % inspect.currentframe().f_code.co_name)

    # initialize new Sgts
     
    new_sgts = Sgts(read=False)
     
    # create output directory
     
    if isinstance(outdir,str):

        if os.path.isfile(outdir):
            ERROR( ("%s file already exists " % outdir) ,exit=True )
            sys.exit()


        if not os.path.isdir(outdir):
            VERBOSE("creating directory %s" % outdir )
            try:
                os.mkdir(outdir)
            except:
                ERROR(("Could not create %s. Exiting") % outdir, exit=True)

        # output file names
         
        out_cgps = outdir+'/vel_cgps.dat'
        out_sgps = outdir+'/vel_sgps.dat'
 
        out_cgps_up = outdir+'/vel_cgps_up.dat'
        out_sgps_up = outdir+'/vel_sgps_up.dat'
         
        # open warning file
         
        fwarning = open(outdir+'/warning.dat','w+')
     
    # start loop on sites
     
    lcode = self.lcode()
     
    for site in lcode:
        VERBOSE("Processing %s" % site)
         
        # check whether there is at least one year of data
         
        if ( self.__dict__[site].data[-1,0] - self.__dict__[site].data[0,0] ) < 1.0 :
            VERBOSE("Less than one year of data for site: %s. Skipping this site." % site)
            if isinstance(outdir,str):
                fwarning.write("-- Less than one year of data for site: %s\n" % site)
            continue

        # check whether there are at least three data
         
        if ( self.__dict__[site].data.shape[0] ) < 3 :
            VERBOSE("Less than 3 data for site: %s" % site)
            if isinstance(outdir,str):
                fwarning.write("-- Less than 3 data for site: %s\n" % site)

        detrended = self.__dict__[site].detrend_median( auto=True )

        if detrended.velocity is None:
            WARNING("Problem in detrend_median for site : %s " % site )
            continue

        if isinstance(outdir, str):
            detrended.save_velocity(out_cgps,verbose=verbose)
            detrended.save_velocity(out_cgps_up,verbose=verbose,up=True)
            new_sgts.append( detrended )
     
    return( new_sgts )       
