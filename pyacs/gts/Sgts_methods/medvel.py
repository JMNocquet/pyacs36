###################################################################
def medvel(self,outdir=None,verbose=False):
###################################################################
    """
     
    Automatic velocity estimates using median estimator.
    The code is adapted from the MIDAS approach (Blewitt et al., 2016).
     
    medvel fills the velocity attribute of every Gts from the current Sgts instance.
     
    returns the modified Sgts instance
    Optionally, if outdir option is provided, writes the results in outdir
     
    :param: outdir: output directory, default None
    :param: verbose: boolean, verbose mode
    :param: warning: output warning file

    :reference: Blewitt, G., Kreemer, C., Hammond, W. C., & Gazeaux, J. (2016). MIDAS robust trend estimator for accurate GPS station velocities without step detection. Journal of Geophysical Research: Solid Earth, 121(3), 2054-2068.
    """

    # import    
    from pyacs.gts.Sgts import Sgts
    import sys , os

     
    # initialize new Sgts
     
    new_sgts = Sgts(read=False)
     
    # create output directory
     
    if outdir is not None:
         
        if os.path.isdir(outdir):
            print('!!! ',outdir,' directory already exists')
            sys.exit()
         
        elif os.path.isfile(outdir):
            print('!!! ',outdir,' file already exists')
            sys.exit()
         
        else:
            if verbose:
                print('-- creating ' , outdir )
             
            os.mkdir(outdir)
         
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

         
        if verbose:
            print('-- Processing ', site )

        # check whether there is at least one year of data
         
        if ( self.__dict__[site].data[-1,0] - self.__dict__[site].data[0,0] ) < 1.0 :
            if verbose:
                print("-- Less than one year of data for site: %s" % site)
            if outdir is not None:
                fwarning.write("-- Less than one year of data for site: %s\n" % site)
            continue

        # check whether there are at least three data
         
        if ( self.__dict__[site].data.shape[0] ) < 3 :
            if verbose:
                print("-- Less than 3 data for site: %s" % site)
            if outdir is not None:
                fwarning.write("-- Less than 3 data for site: %s\n" % site)
#                continue

        # 

        detrended = self.__dict__[site].detrend_median( auto=True )

        if ( detrended.velocity is None ):
            print( "-- Problem in detrend_median for site : %s " % site )
            continue

        if  ( outdir is not None ):

            detrended.save_velocity(out_cgps,verbose=verbose)
            detrended.save_velocity(out_cgps_up,verbose=verbose,up=True)
            new_sgts.append( detrended )
     
    return( new_sgts )       
