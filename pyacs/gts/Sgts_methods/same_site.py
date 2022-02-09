###################################################################
def same_site(self,dc=10, in_place=True, verbose=False):
###################################################################
    """
     
    Check that all gts in the current Sgts are actually the same site. If a given time series is
    found to be of two separate sites, then a new gts is added to the return Sgts instance.

    param dc: critical distance to decide to split the time series
    param in_place: if True modify current Sgts, False retuen a new Sgts
    param verbose: verbose mode
     
    return: a new Sgts instance
    """
     
    # import
    import numpy as np
    from pyacs.gts.Sgts import Sgts
    from pyacs.gts.Gts  import Gts



    if not in_place:
        new_Sgts = Sgts(read=False)

    # start loop on sites
     
    lcode = self.lcode()
     
    for site in lcode:
         
        if verbose:
            print('-- Processing ', site )

        my_ts = self.__dict__[site].copy()

        if my_ts.data_xyz is not None:
            data=my_ts.data_xyz[:,1:4]
            ddata=np.copy(my_ts.data_xyz[:,1:4])
         
        else:
            # if no data_xyz go to next gts
            print("!!! WARNING: data_xyz attribute required for method same_site and not found gts %s" % (site))
         
        # ensure median calculation
        if np.mod(data.shape[0],2)==0:
            # duplicates the last date
            ddata = np.vstack((ddata,ddata[-1,:]))
         
        median = np.median(ddata,axis=0)
        dist_data= np.sqrt( np.sum( (data-median)**2,axis=1) )
         
        lindex = np.where(dist_data > dc*1.E3 )
         
        # case gts needs to be split
        if len( lindex[0] ) > 0 :
            # create a new code
             
            new_code = my_ts.code[:3]+'_'
            if new_code in self.lcode():

                print("!!! ERROR: try to create a new gts with code %s and it already exists." % (new_code))
                new_code = my_ts.code[:2]+'__'
            if verbose:
                print("-- time series for site %s appears to include different sites because there are coordinates at %d dates %.1lf km from the median position" % ( site, len(lindex) , np.max( ddata )*1.E-3 ) )
                print("-- %s time series will be split into code %s and code %s" % (site,site,new_code) )
             
            # create a new gts
            new_gts = Gts(code=new_code,data_xyz=np.copy(my_ts.data_xyz[lindex]))
            new_gts.xyz2neu(corr=True)

            # remove the line from my_ts                

            my_ts.data_xyz = np.delete( my_ts.data_xyz , lindex , axis=0 )

            my_ts.xyz2neu(corr=True)
             
            # update the ouput
             
            if in_place:
                self.append(new_gts)
            else:
                new_Sgts.append( new_gts )
         

        if in_place:
            self.__dict__[site] = my_ts
        else:
            new_Sgts.append( my_ts )

    if in_place:
        return self
    else:
        return new_Sgts
