
import numpy as np
#from ..Gts import Gts
import pyacs

def print_no_return(str):

    import sys
    sys.stdout.write(" %s" % str)
    sys.stdout.flush()


###########################
def __fmt_date(decyear):
    import pyacs.lib.astrotime
    
    date_cal = pyacs.lib.astrotime.decyear2datetime(decyear).strftime("%Y-%m-%d %H:%M")
    dayno,_ut = pyacs.lib.astrotime.decyear2dayno(decyear)

    return(date_cal + (" doy: %03d"%dayno))
###########################


###########################
def __nargmax(data , n , threshold):
    """
    For a given 1D numpy array, returns the index of the n largest values being larger than threshold
    
    return: a 2D numpy array, first column = index, second column = value, line ordered by decreasing value
    """

    import numpy as np

    D = np.copy( data )

    i = 1
    value = 9.E9
    
    lindex = []
    lvalue = []
    
    Dmin = np.min( D )
    
    while ( ( i <= n ) and ( value >= threshold) ):

        
        arg_largest = np.argmax( D )
        value = D[arg_largest]
        
        
        if value >= threshold:
            lindex.append( arg_largest )
            lvalue.append( value )
            
            i = i + 1
            
            D[arg_largest] = Dmin
        
    return  np.array( [ lindex , lvalue ] ).T 

###############################################################################
def suspect_offsets_mf(self, threshold=3,verbose=True,lcomponent='NE',n_max_offsets=5,in_place=False, debug = False):
###############################################################################
    """
    Tries to find offsets in a time series using a median filter
    """
    
    #########################
    def __print_suspects(D):
        for i in np.arange(D.shape[0]):
            print("%10.6lf %04d %8.1lf " % ( D[i,0] , D[i,1] , D[i,2]))
    
    # number of offsets to be searched simultaneously
    n = n_max_offsets
    
    # window lengths (days) for search
    lwindow=[3,5,11,13,15,17,19,21,23,25,27,29,31,33]
    
    suspects_all_date = None
    
    for wl in lwindow:
        
        tmp_ts = self.median_filter(wl).differentiate()
        
        #self.plot(superimposed=self.median_filter(wl).median_filter(wl))
        
        abs_data = np.abs( tmp_ts.data[:,1:4] ) 
        norm_abs_data = abs_data / np.median( abs_data , axis=0 ) 
   
        #####################     
        if 'N' in lcomponent:
            suspects_N =   __nargmax( norm_abs_data[:,0] , n , threshold)
            
            if debug:
                print("-- North wl=%d "% wl)
            if suspects_N.size > 0 :
                ldate = tmp_ts.data[ suspects_N[:,0].astype(int) , 0 ]
                suspects_N_date = np.insert( suspects_N, 0, ldate , axis=1)
                
                if suspects_all_date is None:
                    suspects_all_date = suspects_N_date
                else:
                    suspects_all_date = np.vstack( ( suspects_all_date , suspects_N_date ) )
                    
                
                if debug:
                    __print_suspects(suspects_N_date)
                
            else:
                if debug:
                    print('  -- no suspect found for component North')

        #####################     
        if 'E' in lcomponent:
            suspects_E = __nargmax( norm_abs_data[:,1] , n , threshold)
            
            if debug:
                print("-- East wl=%d "% wl)
            if suspects_E.size > 0 :
                ldate = tmp_ts.data[ suspects_E[:,0].astype(int) , 0 ]
                suspects_E_date = np.insert( suspects_E, 0, ldate , axis=1)

                if suspects_all_date is None:
                    suspects_all_date = suspects_E_date
                else:
                    suspects_all_date = np.vstack( ( suspects_all_date , suspects_E_date ) )
                
                if debug:
                    __print_suspects(suspects_E_date)
                
            else:
                if debug:
                    print('  -- no suspect found for component East')

        #####################     
        if 'U' in lcomponent:
            suspects_U = __nargmax( norm_abs_data[:,2] , n , threshold)
            
            np.save('test.npy',norm_abs_data[:,2])
            
            if debug:
                print("-- Up wl=%d "% wl)
            if suspects_U.size > 0 :
                ldate = tmp_ts.data[ suspects_U[:,0].astype(int) , 0 ]
                suspects_U_date = np.insert( suspects_U, 0, ldate , axis=1)

                if suspects_all_date is None:
                    suspects_all_date = suspects_U_date
                else:
                    suspects_all_date = np.vstack( ( suspects_all_date , suspects_U_date ) )

                if debug:
                    __print_suspects(suspects_U_date)
                
            else:
                if debug:
                    print('  -- no suspect found for component Up')

    ######################
    if verbose:
        print('-- merging the information and selecting most obvious suspect offsets')
    
    if suspects_all_date is not None:
        suspects_all_date=suspects_all_date[np.argsort(-suspects_all_date[:,2])]        
        
        lindex = []
    
        for i in np.arange( suspects_all_date.shape[0] ):
            
            if int(suspects_all_date[i,1]) not in lindex:
                lindex.append( int(suspects_all_date[i,1]) )
                if len( lindex ) > n :
                    break
                
        loffsets_dates = tmp_ts.data[ lindex , 0]
    
    else:
        loffsets_dates = []
        
    if verbose:
        print('-- Suspected offsets (median filter - differentiation) at dates ')
        for date in loffsets_dates:
            print("%10.6lf %s" % (date,__fmt_date(date)) )
    
    ######################

    new_Gts=self.copy()

    if in_place:
        self.offsets_dates = loffsets_dates
        return(self)
        del new_Gts
    else:
        new_Gts=self.copy()
        new_Gts.offsets_dates = loffsets_dates
        return(new_Gts)
  
        

###############################################################################
def test_offset_significance(self,date,conf_level=95,lcomponent='NE',verbose=True,debug=False,mode='local'):
###############################################################################
    """
    test the significance of an offset
    
    :param: date        : date of the offset to be tested
    :param: conf_level  : confidence level in percent (default=95)
    :param: lcomponent  : component to be tested (default='NE')
    :param: mode        : choose among 'local','detrend','detrend_seasonal' to test significance
    :param: verbose     : verbose mode
    :return:            : True if significant, else False
    """
        
    
    if verbose:
        print("-- Testing offset for site %s date %s using mode %s on components %s with confidence level %6.1lf %%" % (self.code,__fmt_date(date),mode,lcomponent,conf_level))


    # F_ratio test
    ##############################################
    def f_ratio(chi_square_1,p1,chi_square_2,p2,n):
        """
        returns result of a F_ratio test
        """
        F=( (chi_square_1-chi_square_2)/(p2-p1) ) / (chi_square_2 / (n-p2) )
        
        from scipy.stats import f
        return(f.cdf(F,p2-p1,n-p2))
        


    
    OK=False
    component={}
    component[1]='North'
    component[2]='East'
    component[3]='Up'
    
    li=[]
    if 'N' in lcomponent:li.append(1)
    if 'E' in lcomponent:li.append(2)
    if 'U' in lcomponent:li.append(3)
    
    H={}
    score={}


    ## TEST LOCAL OFFSET USING A SMOOTHING STRATEGY
    
    if mode=='smooth':

        # local approach uses a weighted value of the score obtained for different time windows around the date
        # number of days before and after the tested date
        l_n_days =  [5,7,9,11,13,15,17]
        lwindow_len=[3,3,5,5 , 5, 7, 7]
        # associated weight to compute the statistics
        coeff=[0.5,1,3,4,10,8,6]


        for i in sorted(li):
        
            print_no_return("-- n_samples/probability %5s:" % component[i] ) 
            
            k=-1
            for n in l_n_days:
                k = k+1
                before = np.array([])
                after = np.array([])
                around = np.array([])
                
                try:
                    before=self.extract_ndates_before_date(date,n).smooth(window_len=lwindow_len[k]).data[:,i] - self.extract_ndates_before_date(date,n).data[:,i]
                    after=self.extract_ndates_after_date(date,n).smooth(window_len=lwindow_len[k]).data[:,i] - self.extract_ndates_after_date(date,n).data[:,i]
                    around=self.extract_ndates_around_date(date,n).smooth(window_len=lwindow_len[k]).data[:,i] - self.extract_ndates_around_date(date,n).data[:,i]
                except:
                    if debug:
                        print("!!! Could not test offset significancy for date %lf using n=%d surrounding data" % (date,n) )
                    break
                
                if debug:
                    print("-- Testing significancy for date %lf using n=%d surrounding data" % (date,n) )

                if before.shape[0] * after.shape[0] * around.shape[0] != 0.0:
    
                    before_std=np.std(before) * 1.E3
                    after_std=np.std(after) * 1.E3
                    around_std=np.std(around) * 1.E3
                    
                    
                    chi_square_2 = before.shape[0] * before_std**2 + after.shape[0] * after_std**2
                    chi_square_1 = around.shape[0] * around_std**2
                    
                    H[i,k]=f_ratio(chi_square_1,1,chi_square_2,2,around.shape[0])*100.0
                    print_no_return(" %02d %5.2lf%% " % (n,H[i,k]))
            
            print("")        
        
            
        
        for i in sorted(li):

            summ=0.0
            score[i]=0.0

            if H != {}:
            
                k=-1
                for n in l_n_days:
                    k=k+1
                    try:
                        score[i]+=H[i,k]*coeff[k]
                        summ=summ+coeff[k]
                    except:
                        if debug:
                            print("-- No data for n=%d " % n)
                        pass
                    
                score[i]=score[i]/summ

            if verbose:
                print("-- Probability of %s offset for site %s date %s component %5s using model: %6.1lf %%" \
                      % (mode.upper(),self.code,__fmt_date(date),component[i],score[i]))


    ## TEST LOCAL OFFSET
    
    if mode=='local':

        # local approach uses a weighted value of the score obtained for different time windows around the date
        # number of days before and after the tested date
        l_n_days = range(2,10)
        # associated weight to compute the statistics
        coeff=[0.5,1,1,4,10,8,6]


        for i in sorted(li):
        
            print_no_return("-- n_samples/probability %5s:" % component[i] ) 

            for n in l_n_days:
                
                before = np.array([])
                after = np.array([])
                around = np.array([])
                
                try:
                    before=self.extract_ndates_before_date(date,n).data[:,i]
                    after=self.extract_ndates_after_date(date,n).data[:,i]
                    around=self.extract_ndates_around_date(date,n).data[:,i]
                except:
                    if debug:
                        print("!!! Could not test offset significancy for date %lf using n=%d surrounding data" % (date,n) )
                    break
                
                if debug:
                    print("-- Testing significancy for date %lf using n=%d surrounding data" % (date,n) )

                if before.shape[0] * after.shape[0] * around.shape[0] != 0.0:
    
                    before_std=np.std(before) * 1.E3
                    after_std=np.std(after) * 1.E3
                    around_std=np.std(around) * 1.E3
                    
                    
                    chi_square_2 = before.shape[0] * before_std**2 + after.shape[0] * after_std**2
                    chi_square_1 = around.shape[0] * around_std**2
                    
                    H[i,n]=f_ratio(chi_square_1,1,chi_square_2,2,around.shape[0])*100.0
                    print_no_return(" %02d %5.2lf%% " % (n,H[i,n]))
            
            print("")        
        
        for i in sorted(li):

            summ=0.0
            score[i]=0.0
            
            for n in range(2,7):
                
                try:
                    score[i]+=H[i,n]*coeff[n]
                    summ=summ+coeff[n]
                except:
                    if debug:
                        print("-- No data for n=%d " % n )
                    pass
            score[i]=score[i]/summ

            if verbose:
                print("-- Probability of %s offset for site %s date %s component %5s using model: %6.1lf %%" \
                      % (mode.upper(),self.code,__fmt_date(date),component[i],score[i]))

    ## TEST GLOBAL OFFSET DETREND
    if mode=='detrend':
        for i in sorted(li):
            
            tmp_ts = self.copy()
            tmp_ts.offsets_dates = []
            #if verbose:
            #    if i==1:print "  => Testing component: North"
            #    if i==2:print "  => Testing component: East"
            #    if i==3:print "  => Testing component: Up"
            
            
            chi_square_1=np.std(tmp_ts.detrend().data[:,i])**2 * tmp_ts.data.shape[0]
            chi_square_2=np.std(tmp_ts.add_offsets_dates([date]).detrend().data[:,i])**2 * tmp_ts.data.shape[0]

            n=tmp_ts.data.shape[0]

            score[i]=f_ratio(chi_square_1,2,chi_square_2,3,tmp_ts.data.shape[0])*100.0

            if verbose:
                print("-- Probability of %s offset for site %s date %s component %5s using model: %6.1lf %%" \
                      % (mode.upper(),self.code,__fmt_date(date),component[i],score[i]))

            
    ## TEST GLOBAL OFFSET DETREND SEASONAL
    if mode=='detrend_seasonal':
        for i in sorted(li):
        
            #if verbose:
            #    if i==1:print "  => Testing component: North"
            #    if i==2:print "  => Testing component: East"
            #    if i==3:print "  => Testing component: Up"
            
            chi_square_1=np.sum(self.detrend().data[:,i]**2)
            chi_square_2=np.sum(self.add_offsets_dates([date]).detrend_seasonal().data[:,i]**2)

            n=self.data.shape[0]

            score[i]=f_ratio(chi_square_1,1,chi_square_2,3,2*n)*100.0

            if verbose:
                print("-- Probability of %s offset for site %s date %s component %5s using model: %6.1lf %%" \
                      % (mode.upper(),self.code,__fmt_date(date),component[i],score[i]))

    
    score_max=np.max(np.array(score.values()))
    if score_max>conf_level:OK=True
    
    if verbose:
        if OK:
            print("-- Offset IS significant at the %5.2lf%% confidence level (threshold=%5.2lf%%) " % (score_max,conf_level))
        else:
            print("-- Offset NOT significant %5.2lf%% confidence level (threshold=%5.2lf%%) " % (score_max,conf_level))
    return(OK)

###########################################################################
def local_offset_robust(self,date,n,verbose=False,debug=False):
    """
    estimate a local offset (no velocity) with a robust method
    :param date: date in decimal year
    :param n: number of dates before and after the dates used in the estimation
    :return : a 1D numpy array with [date, north, east, up, s_north, s_east, s_up]
    :note: the offsets are estimated using the difference between the median position before and after the earthquake using i days
    for all i <=n. Then the median of the estimates is returned. 
    """

    tmp_gts = self.extract_ndates_around_date(date,n)
    array_of_offsets_estimates = None
    
    for i in np.arange(1,n):
        #print tmp_gts.extract_ndates_before_date(date,i).data[:,1:4].shape
        #print tmp_gts.extract_ndates_after_date(date,i).data[:,1:4].shape
        
        before = np.median( tmp_gts.extract_ndates_before_date(date,i).data[:,1:4].reshape(-1,3) , axis=0 ) 
        after  = np.median( tmp_gts.extract_ndates_after_date(date,i).data[:,1:4].reshape(-1,3) , axis=0 ) 
        
        offsets_value = after - before
        #print i,offsets_value*1.E3
        
        if array_of_offsets_estimates is None:
            array_of_offsets_estimates = offsets_value.reshape(-1,3)
        else:
            array_of_offsets_estimates = np.vstack( ( array_of_offsets_estimates , offsets_value ) )
    
    if debug:
        print('all offsets value')
        print(array_of_offsets_estimates*1.E3)
    
    final_values = np.median( array_of_offsets_estimates , axis=0 )

    if verbose:    
        print( "-- combined robust offset estimates (mm): %.1lf " %  final_values * 1.E3 )
    
    return np.array( [date , final_values[0], final_values[1], final_values[2], 0.0, 0.0, 0.0 ] )
       
###########################################################################
def apply_offsets(self,np_offset,opposite=False,in_place=False,verbose=False):
###########################################################################
    """
    Applies given offsets to a times series
    np_offset is a 1D np.array with lines [dates,north,east,up]
    if in_place = True then replace the current time series
    
    :param np_offset: 1D or 3D numpy array or list or list of list with column offset_dates, north, east, up, s_north, s_east, s_up
    :param opposite: boolean, if True apply the oppsite of provided offsets
    :param in_place: if True, will make change in place, if False, return s a new time series
    :param verbose: boolean, verbose mode

        
    """

    # import
    import copy
    import pyacs
    
    # copy
    new_Gts=copy.deepcopy(self)

    # check np_offset

    if not isinstance(np_offset,np.ndarray):
        # np_offset was provided as a list or a list of list
        np_offset=np.array(__ensure_list_of_list(np_offset))
        if ( np_offset.ndim != 2 ) :
            raise TypeError("!!! ERROR: Could not understand offset parameter type:" % np_offset)
            return(self)
        else:
            if ( np_offset.shape[1] < 4 ):
                raise TypeError("!!! ERROR: Could not understand offset parameter type:" % np_offset)
                return(self)
                

    # opposite case
    
    if opposite:
        np_offset[:,1:4] = -np_offset[:,1:4] 
    
        
    new_data=self.data.copy()
    cumulated_offset=np.zeros(new_data.shape)
    
    for i in np.arange(np_offset.shape[0]):
    
        date=np_offset[i,0]
        
        if verbose:
            print("-- adding %10.5lf %10.2lf %10.2lf %10.2lf (mm)" % (date,np_offset[i,1]*1000.0,np_offset[i,2]*1000.0,np_offset[i,3]*1000.0))
        
        lindex=np.where(new_data[:,0] > date)
        
        cumulated_offset[:,1][lindex]=cumulated_offset[:,1][lindex]+np_offset[i,1]
        cumulated_offset[:,2][lindex]=cumulated_offset[:,2][lindex]+np_offset[i,2]
        cumulated_offset[:,3][lindex]=cumulated_offset[:,3][lindex]+np_offset[i,3]
        
        cnew_data=new_data + cumulated_offset
        
        new_Gts.data=cnew_data.copy()

    if in_place:
        self.data=new_Gts.data.copy()
        return(self)
    else:
        return(new_Gts)
        
###############################################################################
def find_offsets_t_scan( self, threshold=0.8, window=250, in_place=False, lcomponent='NE' , verbose=True, debug=True):
###############################################################################
    
    # working Gts
    
    tmp_ts = self.copy()
    
    # import
    
    from pyacs.gts.lib import step_detect

    lindex_step_north = []
    lindex_step_east  = []
    lindex_step_up    = []

    #####################     
    if 'N' in lcomponent:
        
        # t statistics
        t_stat  = step_detect.t_scan(tmp_ts.data[:,1], window=window)
        t_stat_max = t_stat.max()
        t_stat /= np.abs(t_stat).max()

        # threshold
        
        lindex_step_north = step_detect.find_steps(np.abs(t_stat), threshold)

        if debug:
            print("-- North suspected: %d significance: %10.1lf" % (len(lindex_step_north),t_stat_max))
            
    #####################     
    if 'E' in lcomponent:
        
        # t statistics
        t_stat  = step_detect.t_scan(tmp_ts.data[:,2], window=window)
        t_stat_max = t_stat.max()
        t_stat /= np.abs(t_stat).max()

        # threshold
        
        lindex_step_east = step_detect.find_steps(np.abs(t_stat), threshold)
        
        if debug:
            print("-- East  suspected: %d significance: %10.1lf" % (len(lindex_step_east),t_stat_max))
    
    #####################     
    if 'U' in lcomponent:
        
        # t statistics
        t_stat  = step_detect.t_scan(tmp_ts.data[:,3], window=window)
        t_stat_max = t_stat.max()
        t_stat /= np.abs(t_stat).max()

        # threshold
        
        lindex_step_up = step_detect.find_steps(np.abs(t_stat), threshold)
        
        if debug:
            print("-- Up    suspected: %d significance: %10.1lf" % (len(lindex_step_up),t_stat_max ) )
        
    
    ######################

    
    if verbose:
        print('-- merging the information and selecting most obvious suspect offsets')

    suspects_all_date = lindex_step_north + lindex_step_east + lindex_step_up
#     print suspects_all_date
#     print suspects_all_date.sort()
#     print set(suspects_all_date.sort())
#     print list(set(suspects_all_date.sort()))
     
    lindex = list(set(sorted(suspects_all_date)))
    
    if lindex != []:

        lindex_1 = (np.array(lindex) - 1).tolist()                
        loffsets_dates = ( tmp_ts.data[ lindex , 0] + tmp_ts.data[ lindex_1 , 0] ) /2.
    
    else:
        loffsets_dates = []
        
    if verbose:
        print("-- Suspected offsets (t-statistics) at dates")
        for date in loffsets_dates:
            print("%10.6lf %s" % (date,__fmt_date(date)) )

    
    ######################

    new_Gts=self.copy()

    if in_place:
        self.offsets_dates = loffsets_dates
        return(self)
        del new_Gts
    else:
        new_Gts=self.copy()
        new_Gts.offsets_dates = loffsets_dates
        return(new_Gts)


###########################################################################

 
###########################################################################
###########################################################################
# OLD METHODS
###########################################################################
###########################################################################
    

###########################################################################
def get_suspected_dates(diff_data,threshold,lcomponent='NEU',verbose=False):
###########################################################################
    """
    get the list of the largest values; these are suspected offsets
    """
    
    ldates=[]
    [median_north,median_east,median_up,sn,se,su]=np.median(np.abs(diff_data),axis=0)[1:7]
    
    if verbose:
        print("-- median of absolute differentiated time series N %10.3lf    E %10.3lf    U %10.3lf    mm " % (median_north*1.E3,median_east*1.E3,median_up*1.E3))
    
    if 'N' in lcomponent:
        lindex_north=np.where(np.abs(diff_data[:,1]) > threshold*median_north)[0].tolist()

        if verbose:
            print("-- found %d potential offsets from day-to-day difference for component North" % len(lindex_north) ) 

        if lindex_north != []:lindex_north=(np.array(check_suspected_offsets(lindex_north,verbose=verbose)),)
        ldates.append(diff_data[lindex_north][:,0].tolist())
        
        if verbose:
            print("-- adding %d potential offsets for component North" % len(lindex_north) )
        
    if 'E' in lcomponent:
        
        lindex_east=np.where(np.abs(diff_data[:,2]) > threshold*median_east)[0].tolist()

        if verbose:
            print("-- found %d potential offsets from day-to-day difference for component East" % len(lindex_east) ) 

        if lindex_east != []:lindex_east=(np.array(check_suspected_offsets(lindex_east,verbose=verbose)),)
        ldates.append(diff_data[lindex_east][:,0].tolist())

        if verbose:
            print("-- adding %d potential offsets for component East" % len(lindex_east) )


    if 'U' in lcomponent:
        lindex_up=np.where(np.abs(diff_data[:,3]) > threshold*median_up)[0].tolist()
        
        print("-- found %d potential offsets from day-to-day difference for component Up" % len(lindex_up) ) 
        
        
        if lindex_up != []:lindex_up=(np.array(check_suspected_offsets(lindex_up,verbose=verbose)),)

        ldates.append(diff_data[lindex_up][:,0].tolist())

        if verbose:
            print("-- adding %d potential offsets for component Up" % len(lindex_up) )

    lldates=sorted(list(set(ldates[0])))
    return(lldates)


###########################################################################
def check_suspected_offsets(lindex , verbose=False):
###########################################################################
    """
    Check that the list of suspected index does not contain two successive values.
    In this case, it is certainly an isolated outlier and the suspected offsets are removed from the list.
    """
    
    if verbose:
        print('-- checking successive index as indicatoin for outliers rather than offset')
    
    diff_lindex=np.diff(np.array(lindex))
    lbad_index=[]
    for i in np.arange(len(diff_lindex)):
        if diff_lindex[i]==1:
            lbad_index+=[i,i+1]

    if verbose:
        print("-- found %d successive index" % len(lbad_index) )

    
    lgood_index=np.delete(lindex,lbad_index)
    return(lgood_index)

###############################################################################
def __ensure_list_of_list(ll):
###############################################################################
    """
    Ensures ll is a list of lists
    e.g. [a,b] returns [[a,b]],[[a,b]] returns [[a,b]]
    """

    if not isinstance(ll[0],list):
        return([ll])
    else:
        return(ll)


###############################################################################
def find_offsets(self,threshold=3,n_max_offsets=9,conf_level=95,lcomponent='NE',verbose=True,in_place=False):
###############################################################################
    """
    A simple empirical procedure to find offsets.
    
    :param threshold:    threshold value for offset preliminary detection
    :param n_max_offset: maximum number of offsets to be detected simultaneously
    :param conf_level:   confidence level for a suspected offset to be accepted
    :param lcomponent: components used for offset detection

    :return: a new Gts instance with offsets_dates and outliers now populated
    """
    
    # working ts
    
    tts=self.copy()
    
    # although we do not search them, the algorithm might find some outliers
    
    loutliers_dates=[]

    ###########################
    def print_dates_4digits(L):
        if isinstance(L,float):L=[L]
        lreturn=[]
        for date in sorted(L):lreturn.append("%9.4lf" % date)
        return lreturn
    ###########################
    
    
    gross_threshold=10

    if verbose:
        print("********************************************************************")
        print("-- %s STEP #1: trying to identify large offsets" % self.code)
        print("********************************************************************")

    tmp_gts_gross=tts.suspect_offsets(threshold=gross_threshold,verbose=True,lcomponent=lcomponent,n_max_offsets=n_max_offsets)
    
    #str_offsets_dates=" ".join(print_dates_4digits(tmp_gts_gross.offsets_dates))
    
    if verbose:
        print("-- %1d suspected large offsets found" % (len(tmp_gts_gross.offsets_dates)))
        for odate in tmp_gts_gross.offsets_dates:
            print("-- %s %12.8lf" % (__fmt_date(odate),odate))
#        print_no_return("-- ")
#        print "\n -- ".join(map(__fmt_date,tmp_gts_gross.offsets_dates))
        
    significant_offsets=[]
    
    for offset_date in tmp_gts_gross.offsets_dates:
        
        if tts.test_offset_significance(offset_date,conf_level=conf_level,lcomponent=lcomponent,verbose=verbose,mode='local') \
           and tts.test_offset_significance(offset_date,conf_level=conf_level,lcomponent=lcomponent,verbose=verbose,mode='detrend'):
            significant_offsets.append(offset_date)
        else:
            # there is potentialy a large outlier
            if verbose: print("=> Since it is not an offset, date %10.4lf is potentially an outlier" % offset_date)
            tmp_gts_gross.find_outlier_around_date(offset_date,conf_level=conf_level,n=3, lcomponent=lcomponent,verbose=verbose)
    
    if verbose:print("=> Large offset search : %1d offsets confirmed, %1d were actually outliers" % (len(significant_offsets),len(tmp_gts_gross.outliers)))
    
        
    tts.offsets_dates+=significant_offsets
    if tmp_gts_gross.outliers != []:
        lindex_outliers=tmp_gts_gross.outliers
        new_outliers_dates=tmp_gts_gross.data[lindex_outliers,0].tolist()
        loutliers_dates=new_outliers_dates
    
    # find gross outliers
    loutliers_gross=tmp_gts_gross.find_outliers_by_residuals(threshold=gross_threshold).outliers
    loutliers_dates+=self.data[loutliers_gross,0].tolist()

    if verbose:print("=> Outliers search : %03d outliers found" % (len(loutliers_dates)))

    
    tmp_gts_gross.plot()
    
    # intermediate offset
    if verbose:
        print("********************************************************************")
        print("-- %s STEP #2: trying to identify intermediate size offsets" % self.code)
        print("********************************************************************")
    
    intermediate_threshold=0.5*(10.+threshold)

    # get estimates from previously found offsets
    
    previous_offsets_values=tmp_gts_gross.remove_outliers().detrend().offsets_values
    
    # work on the already corrected times series
    tmp_gts_intermediate=tmp_gts_gross.remove_outliers().apply_offsets(previous_offsets_values).suspect_offsets(threshold=intermediate_threshold,verbose=True,lcomponent='NE',n_max_offsets=10)
    str_offsets_dates=" ".join(print_dates_4digits(tmp_gts_intermediate.offsets_dates))
    if verbose:
        print("=> %1d potential intermediate offsets found at %s" % (len(tmp_gts_intermediate.offsets_dates),str_offsets_dates))

    significant_offsets=[]

    lpotential_offsets_dates=tmp_gts_intermediate.offsets_dates

    tmp_gts_intermediate.offsets_dates=[]

    
    for offset_date in lpotential_offsets_dates:
        if tmp_gts_intermediate.test_offset_significance(offset_date,conf_level=conf_level,lcomponent=lcomponent,verbose=verbose,mode='local') \
           and tmp_gts_intermediate.test_offset_significance(offset_date,conf_level=conf_level,lcomponent=lcomponent,verbose=verbose,mode='detrend'):
            significant_offsets.append(offset_date)
        else:
            # there is potentialy a large outlier
            if verbose: print("=> Since it is not an offset, date %10.4lf is potentially an outlier" % offset_date)
            tmp_gts_intermediate.find_outlier_around_date(offset_date,conf_level=conf_level,n=3, lcomponent=lcomponent,verbose=verbose)
    
    
    if verbose:print("=> Intermediate size offsets search : %1d offsets confirmed, %1d were actually outliers" % (len(significant_offsets),len(tmp_gts_intermediate.outliers)))
            
    tts.offsets_dates+=significant_offsets

    if tmp_gts_intermediate.outliers != []:
        lindex_outliers=tmp_gts_intermediate.outliers
        new_outliers_dates=tmp_gts_intermediate.data[lindex_outliers,0].tolist()
        loutliers_dates=loutliers_dates+new_outliers_dates

    # find intermediate outliers
    loutliers_intermediate=tmp_gts_intermediate.find_outliers_by_residuals(threshold=intermediate_threshold).outliers
    loutliers_dates+=self.data[loutliers_intermediate,0].tolist()

    if verbose:print("=> Outliers search : %03d outliers found" % (len(loutliers_dates)))


    # final search

    if verbose:
        print("********************************************************************")
        print("-- %s STEP #3: trying to identify small offsets" % self.code)
        print("********************************************************************")

    # get estimates from previously found offsets
    
    previous_offsets_values=tmp_gts_gross.remove_outliers().detrend().offsets_values

    tmp_gts_final=tmp_gts_intermediate.remove_outliers().apply_offsets(previous_offsets_values).suspect_offsets(threshold=threshold,verbose=True,lcomponent='NE',n_max_offsets=10)

    str_offsets_dates=" ".join(print_dates_4digits(tmp_gts_final.offsets_dates))
    if verbose:print("=> %1d potential subtle offsets found at %s" % (len(tmp_gts_final.offsets_dates),str_offsets_dates))

    significant_offsets=[]

    lpotential_offsets_dates=tmp_gts_final.offsets_dates

    tmp_gts_final.offsets_dates=[]

    
    significant_offsets=[]
    for offset_date in lpotential_offsets_dates:
        if tmp_gts_intermediate.test_offset_significance(offset_date,conf_level=conf_level,lcomponent=lcomponent,verbose=verbose,mode='local') \
           and tmp_gts_intermediate.test_offset_significance(offset_date,conf_level=conf_level,lcomponent=lcomponent,verbose=verbose,mode='detrend_seasonal'):
            significant_offsets.append(offset_date)
        else:
            # there is potentialy a large outlier
            tmp_gts_final.find_outlier_around_date(offset_date,conf_level=conf_level,n=3, lcomponent=lcomponent,verbose=verbose)
            

    if verbose:print("=> Subtle offset search: %1d offsets confirmed, %1d were actually outliers" % (len(significant_offsets),len(tmp_gts_final.outliers)))

        
    tts.offsets_dates+=significant_offsets

    # ensure offsets are listed once

    offsets_dates=sorted(list(set(tts.offsets_dates)))
    
    tts.offsets_dates=offsets_dates
    
    # populate outliers

    if tmp_gts_final.outliers != []:
        lindex_outliers=tmp_gts_final.outliers
        new_outliers_dates=tmp_gts_final.data[lindex_outliers,0].tolist()
        loutliers_dates=loutliers_dates+new_outliers_dates

    # final outliers
    loutliers_final=tmp_gts_final.find_outliers_by_residuals(threshold=threshold).outliers
    loutliers_dates+=self.data[loutliers_final,0].tolist()

    if verbose:print("=> Outliers search : %03d outliers found" % (len(loutliers_dates)))
    
    if loutliers_dates != []:
        from ..Gts import get_index_from_dates
        returned_index=get_index_from_dates(loutliers_dates,tts.data,tol=0.25)
        loutliers=sorted(list(set(returned_index)))
        #tts.outliers=loutliers
    else:
        loutliers=[]
    
    if verbose:print("**********************************************************************************")
    if verbose:print("=> Final results of offset search: %02d offsets found; additionally %02d outliers were flagged" % (len(offsets_dates),len(loutliers_dates)))
    if verbose:print("**********************************************************************************")

    new_Gts=self.copy()
    new_Gts.offsets_dates=offsets_dates
    new_Gts.outliers=loutliers

    if in_place:
        self.offsets_dates=new_Gts.offsets_dates
        self.outliers=new_Gts.outliers
        return(self)
        del new_Gts
    else:
        return(new_Gts)

###############################################################################
def suspect_offsets(self, threshold=3,verbose=True,lcomponent='NE',n_max_offsets=10,in_place=False):
###############################################################################
    """
    Tries to find offsets in a time series
    """
    
    # differentiate, returns a time series with dates between points
    
    dates=self.data[:,0]
    diff_data=np.diff(self.data,n=1,axis=0)
    
    new_dates=dates[:-1]+diff_data[:,0]/2.
    
    diff_data[:,0]=new_dates
    
    
#    [dd,median_north,median_east,median_up,sn,se,su]=np.median(np.abs(diff_data),axis=0)
#    if verbose:print("median_north %8.2lf median_east %8.2lf median_up %8.2lf" %(median_north*1000.0,median_east*1000.0,median_up*1000.0))
    

    
    
    lldates=get_suspected_dates(diff_data,threshold , verbose=verbose)
    

    if len(lldates)>n_max_offsets:
        if verbose:
            '-- Too many offsets detected. Keeping only the largest for analysis and will see later the others...'
        while len(lldates)>n_max_offsets:
            threshold=threshold*1.5
            lldates=get_suspected_dates(diff_data,threshold)
        

    new_Gts=self.copy()

    if in_place:
        self.offsets_dates=new_Gts.offsets_dates
        return(self)
        del new_Gts
    else:
        new_Gts=self.copy()
        new_Gts.offsets_dates=lldates
        return(new_Gts)



    


###############################################################################
def find_time_offsets(self,option=None,ndays=7,th_detection_rms=3,th_detection_offset=3):
###############################################################################
    
    """
     Find the time of suspected offsets by rms time series calculated over ndays
     Then check the time of offsets: if one offset is too small/None then it is removed
     input:

    - ndays: number of positions in the windows. rms time series are calculated over ndays.
    - th_detection_rms: the threshold for detecting the anomalous windows rms(t) > th_detection_rms*median(rms(ts)).
    - th_detection_offset: the threshold for detecting the offsets,
       for each anomalous time windows, differentiate positions
       then test whether it is a suspected offset (differentiated(t) > threshold * median(differentiated))

     output:

         add the time of offsets in to self.offsets
    """
    if not option:option='detrend'
    
    ## clean data before using
    data_saved = self.copy().data
    outliers_saved = self.outliers
    if self.outliers: self.remove_outliers(in_place=True)
    ## calculate rms of windows 
    data_rms = self.rms(ndays)     
    
    ##########################################################
    ## find time of possible offsets by rms: anomalous windows
    ##########################################################
    lindex_offset = []
    for i in (1,2,3):
        
        threshold_rms = th_detection_rms*np.median(data_rms[:,i])
        lindex_offset_i = []
        if max(data_rms[:,i]) > threshold_rms:
            
            index_out_rms = np.argwhere(data_rms[:,i] > threshold_rms)                
            
            for j in index_out_rms:                   
                
                differentiate_ndays = np.hstack(np.fabs(np.diff(self.data[j:(j+ndays),i])))
                index = np.argwhere(differentiate_ndays > th_detection_offset*np.median(differentiate_ndays))
                
                ## offset is in the windows what the length of out_index is 1
                if (len(index) == 1): lindex_offset_i.append(index[0] + j)

            if len(lindex_offset_i) > 0:
                # lindex_offset.append(np.hstack(lindex_offset_i))
                lindex_offset_i = np.hstack(lindex_offset_i)
                for i in np.unique(lindex_offset_i):
                    if len(np.argwhere(lindex_offset_i==i)) >= 2: lindex_offset.append(i) ##1: ok for Taiwan
    ## take the time of offsets from time series
    if len(lindex_offset) > 0:
        self.offsets = list(self.data[np.unique(lindex_offset),0])
    else: self.offsets = []
#        ## choose one in two continuous offsets if the number of position between them < ndays
#        if len(t0) > 1: t0 = pygps_module.check_middle_two_continuous_offets(ndays,t0,self.data[:,0])
#        ### check the offsets by their amplitudes
#        if len(t0)>0: self.offsets = self.delete_small_offsets(t0)
    
    ##########################################################
    ## test the offsets found by F_ratio test
    ##########################################################
    ## re-turn to initial data     
    self.data = data_saved
    self.test_offsets()       
    self.outliers = outliers_saved   


###############################################################################
def delete_small_offsets(self,offsets,del_by_pricise=False):
###############################################################################
    """
        The aim for test_offset modul. 
        Estimate the offsets with clean data.
        Then delete the offsets which their values are so small
        input: list of time offsets
        output: list of time offsets tested
    """
    ### make the threshold of offset
    min_offset_NE = 0.002
    if self.H_conf:
        if self.H_conf.has_key('threshold_smallest_offset'):min_offset_NE = self.H_conf['threshold_smallest_offset']
    min_offset_UP = 3.*min_offset_NE
    self.offsets = offsets
    ## estimate the offsets to check its by their amplitude
    new_Gts = self.make_model(loutlier=self.outliers)
    if (new_Gts is not None) and (new_Gts.offsets_values is not None):
        offsets_values = np.fabs(new_Gts.offsets_values)
        ### remove the very small offsets, incorrect offsets
        choose_index = []
        for i in range(len(self.offsets)):
            if del_by_pricise:
                if (offsets_values[i,1] > min_offset_NE) and (offsets_values[i,1] > offsets_values[i,4]): choose_index.append(i)
                if (offsets_values[i,2] > min_offset_NE) and (offsets_values[i,2] > offsets_values[i,5]): choose_index.append(i)
                if (offsets_values[i,3] > min_offset_UP) and (offsets_values[i,3] > offsets_values[i,6]): choose_index.append(i)
            else:
                if (offsets_values[i,1] > min_offset_NE): choose_index.append(i)
                if (offsets_values[i,2] > min_offset_NE): choose_index.append(i)
                if (offsets_values[i,3] > min_offset_UP): choose_index.append(i)
        if len(choose_index)>0:
            choose_index = np.unique(choose_index)
            offsets = list(np.take(self.offsets, choose_index, axis=0))
        else: offsets = None
    self.offsets = offsets
    self.offsets_values = None
    
    return offsets




###########################################################################################################
### F ratio test
### f not defined pb to be solved
###########################################################################################################
# def test_offsets_significance(self,data,t0):
#         
#     """This method makes a test of significancy of an offset in a time series
#        The time series is assumed to have ONLY ONE offset provided at  
#        provided offset using an F_ration test
#        The chi2 of two estimation are compared:
#        - one without offset
#        - one with one offset
#        The F_ratio is then formed to decide wether the offset is significant or not"""
#        
#     n_unknown1 = 2      #without offset
#     n_unknown2 = 2 + 1  # with offset
# 
#     p1=data.shape[0]-n_unknown1
#     p2=data.shape[0]-n_unknown2
# 
#     alpha=0.05
#     #f_distribution = f.ppf((1-alpha),(p1-p2),p2)
#     f_distribution = f.ppf((1-alpha),p2,(p1-p2))
# 
#     self.data = data
#     
#     ## without offset
#     self.offsets = []
#     new_Gts1 = self.make_model(option='detrend')
#     ## with offset
#     self.offsets = [t0]
#     new_Gts2 = self.make_model(option='detrend')
#     
#     if (new_Gts1 is not None) and (new_Gts2 is not None): t0 = None
#     elif (new_Gts1 is not None) and (new_Gts2 is not None):
#         
#         V1 = new_Gts1.data
#         V2 = new_Gts2.data
#         
#         t0_checked = []
#         for i in range(1,4):
#             K1 = np.dot(np.transpose(V1[:,i]),V1[:,i])
#             K2 = np.dot(np.transpose(V2[:,i]),V2[:,i])
#             F_ratio=((K1-K2)/(p1-p2))/(K2/p2)
#             if F_ratio > f_distribution: t0_checked.append(t0)
#             #print t0, F_ratio, f_distribution
#         if len(t0_checked) > 0: t0 = t0
#         else: t0 = None
#     else: t0 = None
# 
#     return t0
    
    
###########################################################################################################
### Test the offsets
###########################################################################################################
def test_offsets(self,verbose=False,debug=True,window_length=None):
    """
        Test the offset:
         - delete the small offset (1mm for East/North, 2mm for Up)
         - then make a F ratio test
         - re-check to delete the small offsets
    """
    
    ## find the outliers before test
    ouliers_saved = self.outliers
    self.find_outliers_by_residuals()
    
    ### TEST BY VALUE OF OFFSETS: estimate the offsets to check its by their amplitude
    self.offsets = self.delete_small_offsets(self.offsets)
    ### TEST BY F_RATIO
    if self.offsets is not None:
        #if not window_length: window_length = 50 ##4months: ok for synthetic data, 2months: ok for INGR
        if not window_length: window_length = 180 ##4months: ok for synthetic data, 2months: ok for INGR
        ### save the original data 
        data_saved = self.copy().data
        ### clean data to test
        data = self.remove_outliers().data            
        t0 = self.offsets

        list_index  = np.hstack([0,np.searchsorted(data[:,0],t0),len(data[:,0])])
        t0_tested = []
        for i in range(len(t0)):

            if verbose:
                print("  -> Testing offset time %10.3lf" % t0[i])
            
            if window_length:
                ## start index of data for test t0_i
                if (list_index[i+1] - list_index[i] < window_length/2): index_start = list_index[i]
                else: index_start = list_index[i+1] - window_length/2
                ## end index of data for test t0_i
                if (list_index[i+2] - list_index[i+1] < window_length/2): index_end = list_index[i+2]
                else: index_end = list_index[i+1] + window_length/2
            else:
                index_start = list_index[i]
                index_end = list_index[i+2]
            
            ## take data_i for test
            data_i = data[index_start:index_end]

            ## test by F test ratio
            t0_i = self.Ftest_4_offsets(data_i,t0[i])
            
            if t0_i is not None:
                if verbose: 
                    print("    -> Significant. Goes to next test.")
                t0_tested.append(t0[i])
            else:
                if verbose: 
                    print("    -> No nignificant. Goes to next test.")
        
        if len(t0_tested) > 0: self.offsets = list(np.unique(t0_tested))
        else: self.offsets = None
        
        self.offsets_values = None
        ## return self.data
        self.data = data_saved 
    
    ### RE-TEST BY VALUE OF OFFSETS
    if self.offsets is not None:
        ## find the outliers before test
        self.find_outliers_by_residuals()
        ## estimate the offsets to check its by their amplitude
        self.offsets = self.delete_small_offsets(self.offsets,del_by_pricise=True)

    self.outliers = ouliers_saved

###########################################################################
## estimate the local offsets of time series with the trendline y = a  + offsets
###########################################################################
def estimate_local_offset(self,window_length=4, in_place=False):
    """
    Estimate the local offset, just used window_length positions before & window_length positions behind of offset
    output: amplitude of local offsets
    """
    
    import pyacs.gts.lib.gts_estimators
    
    ### take the data
    clean_data = self.remove_outliers().data
    
    local_offsets_values = np.zeros((len(self.offsets_dates),7))
    local_offsets_values[:,0] = self.offsets_dates

    
    noffset = len(self.offsets_dates)
    
    list_index_t0_i  = np.searchsorted(clean_data[:,0],local_offsets_values[:,0])        

    index_local = []
    for i in list_index_t0_i:
        index_local.append(range(i-window_length,i+window_length))
    
    index_local = np.hstack(index_local)
    index_local = np.delete(index_local,np.argwhere(index_local<0))
    index_local = np.delete(index_local,np.argwhere(index_local>=len(clean_data[:,0])))

    data_local = np.take(clean_data,index_local,axis=0)
            
    ndate = len(data_local[:,0])
    for k in range(1,4):            
        A = np.ones((ndate,(1+noffset)), float)
        for i in range(ndate):
            for j in range(noffset):
                if data_local[i,0] <= local_offsets_values[j,0]: A[i,(1+j)] = 0.
                
        X,s_X = pyacs.gts.lib.gts_estimators.least_square(A,data_local[:,k])[0:2]
        if X is not None:
            local_offsets_values[:,k] = X[1:]
            local_offsets_values[:,3+k] = s_X[1:]
        else: local_offsets_values = None

    if in_place:
        self.offsets_values=local_offsets_values
        return(self)
    else:
        new_Gts=self.copy()
        new_Gts.offsets_values=local_offsets_values
        return(new_Gts)
    

    return local_offsets_values

