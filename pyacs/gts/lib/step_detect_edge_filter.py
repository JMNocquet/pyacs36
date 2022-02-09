###############################################################################
def find_offsets_edge_filter( self, threshold=0.6, search_lbda=[ 3,5,7,10,20,50,100,200,300] , delta_day=100, \
                              in_place=False, lcomponent='NE' , \
                              verbose=True, debug=True, log=False, eq_file=None):
###############################################################################
    
    # import
    
    from pyacs.gts.lib.filters.total_variation import edge
    import numpy as np
    import pyacs.lib.astrotime as at
    import pyacs.lib.units
    import datetime
    
    # handle delta_day
    
    if isinstance(delta_day,int):
        ldelta_day = [delta_day]
    else:
        if not isinstance(delta_day,list):
            print("!!! ERROR: delta_day should be a list of integer or an integer")
        else:
            ldelta_day = delta_day
    
    # eq list a gcmt file
    if eq_file is not None:
        EQ=np.genfromtxt(eq_file,dtype=None)
        np_eq_lon=EQ['f0']
        np_eq_lat=EQ['f1']
        np_moment_scalar=EQ['f8']
        np_moment_exponent=EQ['f9']
        np_str_date=EQ['f10']
        
        # handle magnitude
        np_eq_magnitude = np.copy( np_moment_scalar ) *0.0
        i = 0
        for scalar,exponent in zip( np_moment_scalar , np_moment_exponent ):
            moment = scalar * 10.**(exponent - 7)
            np_eq_magnitude[i] = pyacs.lib.units.moment_to_magnitude(moment)
            i= i + 1
        # handle dates
        np_eq_mjd = np.copy( np_moment_scalar ) * 0.
        i=0
        for str_date in np_str_date:
            dyear  = int( str_date[:4] )
            dmonth = int( str_date[4:6] )
            dmday  = int( str_date[6:8] )
            dhour  =  int( str_date[8:10] )
            dminute = int( str_date[10:12] )
            date = datetime.datetime(dyear,dmonth,dmday,dhour,dminute)
            np_eq_mjd[i] = at.datetime2mjd( date )
            i = i + 1

    # initialize the components
    l_idx_component = []
    if 'N' in lcomponent:
        l_idx_component.append(1)
    if 'E' in lcomponent:
        l_idx_component.append(2)
    if 'U' in lcomponent:
        l_idx_component.append(3)
    
    H_idx_component = {1:'N',2:'E',3:'U'}
    
    # working Gts
    new_gts = self.copy()

    
    ##############
    def __get_offset_date_and_value( x, y , threshold ):
    
        # for dates
        data_mjd = at.decyear2mjd( x )
        new_data_mjd = ( data_mjd[0:-1] + data_mjd[1:] ) / 2. 
        
        x_diff = at.mjd2decyear( new_data_mjd )
        
        # for y
        
        diff_y = np.diff( y )
        abs_diff_y = np.fabs( diff_y )
        
        np_idx = np.where( abs_diff_y > threshold )
        
        # get dates
        np_dates = x_diff[ np_idx ]
        # get values
        np_values = diff_y[ np_idx ]
        
        return np_dates,np_values
    ##############
    
    ##############
    def __handle_round_date( decyear ):
        mjd = np.round( at.decyear2mjd( decyear ) , decimals=2 )
        (mday,month, iyear,_ut) = at.mjd2cal(mjd)
        return( mday, month, iyear)
    ##############
        
    
    # detrend

    H_offset={}
    
    for delta_day in ldelta_day:
        
        if verbose:
            print("-- detrend_median with delta_day=%d" % delta_day)
    
        gts_dtm = self.detrend_median( delta_day=delta_day )
    
        if gts_dtm is None:
            gts_dtm = self.detrend( )
        
            
        median_neu = np.median( np.fabs( gts_dtm.differentiate().data[:,1:4] ) , axis=0 )
    
        if verbose:
            print("-- median daily scattering NEU (mm): " , median_neu * 1.E3 )
    
            
        for lbda in search_lbda:
            
            edge_ts = gts_dtm.edge_filter(lbda, verbose=verbose)
                
            for idx in l_idx_component:
                
                np_dates,np_values = __get_offset_date_and_value(edge_ts.data[:,0], edge_ts.data[:,idx], median_neu[idx-1]*threshold)
        
                new_gts.offsets_dates = np_dates
                
                if verbose or log:
                    
                    if verbose:                
                        print("-- suspected offsets for component %s and lbda = %d" % (H_idx_component[idx] , lbda ))
    
                    for date, value in zip(np_dates,np_values):
                        ( mday, month, iyear) = __handle_round_date(date)
                        mjd = at.decyear2mjd( date )
                        sstr_date = ("%4d_%02d_%02d" % (iyear,month,mday))
                        
                        no_EQ = True
                        if eq_file is not None:
                            i=0
                            for mjd_eq in np_eq_mjd:
                                if np.fabs( mjd_eq-mjd ) <=1.:
                                    if verbose:
                                        print("%.6lf %10.3lf (%4d-%02d-%02d) EQ Mw %.1lf %6.2lf %6.2lf " % \
                                              (date,value*1.E3,iyear,month,mday, np_eq_magnitude[i],np_eq_lon[i],np_eq_lat[i]))
                                    
                                    if sstr_date in H_offset.keys():
                                        H_offset[sstr_date].append([idx,value,lbda,np_eq_magnitude[i],delta_day])
                                    else:
                                        H_offset[sstr_date] = [[idx,value,lbda,np_eq_magnitude[i],delta_day]] 
    
                                    no_EQ=False
                                i = i + 1
                        if no_EQ and verbose:            
                            print("%.6lf %10.3lf (%4d-%02d-%02d)" % (date,value*1.E3,iyear,month,mday))
    
                        if no_EQ:            
                            if sstr_date in H_offset.keys():
                                H_offset[sstr_date].append([idx,value,lbda,0,delta_day])
                            else:
                                H_offset[sstr_date] = [[idx,value,lbda,0,delta_day]] 

    if log:
        for sstr_date in H_offset.keys():
            # open log file
            
            log_filename=("%s_Mw%02d.dat" % (sstr_date,H_offset[sstr_date][0][-2]*10))
            try:
                flog=open(log_filename,'a+')
            except:
                from colors import red
                return_str = red( "[PYACS WARNING] Could not create file: %s" % ( log ) ) 
                print( return_str )
            print(H_offset[sstr_date])
            np_offset = np.array(H_offset[sstr_date]).reshape(-1,5)
            # list lbda
            l_lbda = list(set(np_offset[:,2].tolist()))
            # list dalta_day
            l_delta_day = list(set(np_offset[:,4].tolist()))
            # loop lbda
            for lbda in sorted(l_lbda):
                # loop delta_day
                for delta_day in l_delta_day:
                    # empty line to write as a list
                    lw = [self.lon,self.lat,0.,0.,1.,1.,0.,("%s_%03d" % (self.code,lbda)),0]
                    # get results given ldba
                    lindex = np.where( (np.fabs(np_offset[:,2] - lbda) < 1.E-3) & ( np.fabs(np_offset[:,4] - delta_day) < 1.E-3 ) )
                    offset = np_offset[lindex].reshape(-1,5)
                    for i in np.arange( offset.shape[0] ):
                        if offset[i,0]==1.:
                            lw[3] = offset[i,1]*1.E3
                        if offset[i,0]==2.:
                            lw[2] = offset[i,1]*1.E3
                        if offset[i,0]==3.:
                            lw[7] = ("%s %6.2lf " % (lw[7]),offset[i,1]*1.E3 )
                    
                        lw[8] = offset[i,4]
                
                flog.write("%10.4lf %10.4lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf %s %03d\n" % tuple(lw) )
            
            flog.close()

    return new_gts
                
