###################################################################
## easy: a basic analysis of time series and plot
###################################################################

def lazy_pyacs(self,verbose=True,in_place=False):
    """
    A simple automated analysis for lazy analysts
    """
    # initialization
    
    import os
    
    conf_level_offset=95.0 # confidence level for accepting an offset detection
    itermax=10 # max iteration for outliers detection
    threshold_outlier=3.8 # threshold for outlier detection (outliers are flagged if > threshold_outlier * median)
    n_swindow=20 # n samples for the sliding windows
    loutliers_dates=[]
    
    
    # print info
    if verbose:
        print("********************************************************************")
        print("lazy_pyacs analysis of time series %s (%s)" % (self.code,self.file))
        print("********************************************************************")
        print("Parameters used for the analysis:")
        print("    - Confidence level for offset detection              : %4.2lf%%" % conf_level_offset)
        print("    - Threshold for outlier detection                    : %4.2lf" % threshold_outlier)
        print("    - Time window length used for outlier detection      : %02d" % n_swindow)
        print("    - Maximum number of iterations for outlier detection : %02d" % itermax)
        print("********************************************************************")
    
    # find offsets

    if verbose:
        print("********************************************************************")
        print("-- %s  STEP #1: trying to identify offsets" % self.code)
        print("********************************************************************")

    tmp_ts = self.suspect_offsets_mf(1000)        
#        tmp_ts=self.find_offsets(conf_level=95.0,verbose=verbose)
    loutliers_find_offsets=tmp_ts.outliers
#        loutliers_dates_find_offsets=tmp_ts.data[loutliers_find_offsets,0].tolist()
#        loutliers_dates+=loutliers_dates_find_offsets
    loffsets_dates=tmp_ts.offsets_dates



    if verbose:
        print("********************************************************************")
        print("-- %s STEP #1: %02d offsets  found" % (self.code,len(loffsets_dates)))
        print("-- %s STEP #1: %02d outliers found" % (self.code,len(loutliers_find_offsets)))
        print("********************************************************************")
    
    # outliers through sliding windows iteratively

    print "********************************************************************"
    print("-- %s STEP #1: Trying to identify outliers using rms and sliding window strategy" % self.code)
    print "********************************************************************"

    # sliding window
    tmp_ts.outliers=[]
    loutliers_sw=tmp_ts.find_outliers_sliding_window(threshold=5.0).outliers
    loutliers_dates_sw=tmp_ts.data[loutliers_sw,0].tolist()
    loutliers_dates+=loutliers_dates_sw
    
    if verbose:        
        print("-- %s STEP #4: %02d additional outliers found using sliding window" % (self.code,len(loutliers_dates_sw)))
    
    # residuals
    tmp_ts.outliers=loutliers_sw
    tmp_ts=tmp_ts.remove_outliers()
    loutliers_res=tmp_ts.find_outliers_by_residuals(threshold=5.0).outliers
    loutliers_dates_res=tmp_ts.data[loutliers_res,0].tolist()
    loutliers_dates+=loutliers_dates_res

    if verbose:        
        print("-- %s STEP #5: %02d additional outliers found using residuals" % (self.code,len(loutliers_dates_res)))

    

    new_Gts=self.copy()
    new_Gts.outliers=get_index_from_dates(loutliers_dates,self.data,tol=0.25)
    new_Gts.offsets_dates=loffsets_dates
    
    if not os.path.isdir('../glred'):
        try:
            os.mkdir('../glred')
        except:
            print '!!! Could not create ../glred'
        
    eq_rename='../glred/'+new_Gts.code+'.eq_rename'
    if verbose:
        print("********************************************************************")
        print("-- Final results of lazy_pyacs analysis: %02d offsets  detected ; %03d outliers detected" % \
              (len(new_Gts.offsets_dates),len(new_Gts.outliers)))
        print("-- Saving outliers and offsets_dates to %s " % eq_rename)
        print "********************************************************************"

    #new_Gts.save_eq_rename(eq_rename,verbose=False, excluded_periods=None)


    if in_place:
        self.outliers=new_Gts.outliers
        self.offsets_dates=new_Gts.offsets_dates
        return(self)
        del new_Gts
    else:
        return(new_Gts)

###############################################################################
def lazy(self , conf_level_subtle=0 , itermax=10, components='NE' , in_place=False, verbose=True):
###############################################################################

    ###########################
    def __fmt_date(decyear):
        import pyacs
        
        date_cal = pyacs.decyear2datetime(decyear).strftime("%Y-%m-%d %H:%M")
        dayno,ut = pyacs.decyear2dayno(decyear)
    
        return(date_cal + (" doy: %03d"%dayno))

    ###########################################################
    def __remove_already_tested_dates(dates,ref_dates, tolerance=0.5):
        
        import numpy as np
        
        lgood_dates = []
        np_ref_dates = np.array(ref_dates)
        
        one_day_dec_year = 1. / 365.25
        
        tol_year = tolerance * one_day_dec_year
        
        for date in dates:
            
            if np.min( np.abs( np_ref_dates - date) ) > tol_year:
                lgood_dates.append( date )
            
        return( lgood_dates )
    
    ###########################################################
    # STEP 1 FIND GROSS OFFSET
    ###########################################################
     
    # days used to estimate local offset
    ndays = 7
    
    # threshold
    
    threshold=200.
    
    loffsets = []
    lnot_signicant_offset = []
    tmp_gts = self.detrend_seasonal_median()
    
    if verbose:
        print '-- STEP #1 GROSS OFFSET SEARCH USING MEDIAN_FILTER AND DIFFERENTIATION WITH THRESHOLD=',threshold

    # BEGIN LOOP
    OK=True
    iteration = 0
    while OK:

        iteration = iteration + 1
        print '-- Current wrms is      = ', tmp_gts.detrend_seasonal().wrms() * 1.E3  
        print '-- Current offset list  = ', loffsets  
        print '-- Running suspect_offsets_mf with threshold = ', threshold , ' iteration #', iteration 

        tmp_gts = tmp_gts.suspect_offsets_mf(threshold,n_max_offsets=10,lcomponent=components,verbose=verbose,debug=False)
#            tmp_gts = tmp_gts.find_offsets_t_scan(threshold=threshold,window=50,lcomponent=components,verbose=verbose,debug=False)


        # Eliminate already found offsets
        print '-- Cheking whether these suspected offsets are new'
        
        if loffsets != []:
            tmp_gts.offsets_dates = __remove_already_tested_dates(tmp_gts.offsets_dates,loffsets)

        print '-- There are ',len(tmp_gts.offsets_dates),' potential offsets to be tested'

        # IF NO OFFSET SUSPECTED GO TO NEXT STEP             
        if tmp_gts.offsets_dates == []:
            if verbose:
                print '-- No more offset suspected for threshold = ',threshold,'. Going to next step.'
            OK=False
            break

        
        # TEST SUSPECTED OFFSETS
        if verbose:
            print '-- found ',len( tmp_gts.offsets_dates ),' potential offsets'
            
        print '-- Testing significance of potential offsets'
        
        ldate = tmp_gts.offsets_dates

        tmp_gts.offsets_dates = []
        
        for date in ldate:
            
            print '-- Testing potential offset at date: ',date,__fmt_date(date)
            
            SIGNIFICANT=tmp_gts.test_offset_significance(date, lcomponent=components,mode='local')
             
            if SIGNIFICANT:
                print '-- offset at date ',date, __fmt_date(date),' is significant'
                loffsets.append( date )
                offset_values = tmp_gts.local_offset_robust(date , ndays ,verbose=verbose)
                print '-- correcting offset at date ',date, __fmt_date(date),' with values ', offset_values[1:4]*1.E3, 'mm'
                tmp_gts = tmp_gts.apply_offsets(offset_values.reshape(-1,7),opposite=True)
                
            else:
                print '-- Offset at date ',date,__fmt_date(date),' not significant'
                lnot_signicant_offset.append(date)

        # test whether the maximum iteration has been reached
        if iteration >= itermax:
            print '-- Maximum number of iteration reached: ', itermax
            OK=False
    # END LOOP

        
    print '-- Kept offset_dates ',loffsets
    
    
    
    # Find Outliers
    tmp_gts.offsets_dates = loffsets
    tmp_gts = tmp_gts.detrend_seasonal().find_outliers_vondrak(threshold=5., fc=2. , verbose=True, component=components)
    
    # Final Global test on offsets

    loutliers = tmp_gts.outliers
    ldate = tmp_gts.offsets_dates

    lfinal_offset_dates = []
    
    tmp_gts = self.copy()
    
    tmp_gts.offsets_dates = []
    tmp_gts.outliers = loutliers
    tmp_gts = tmp_gts.remove_outliers() 

    for date in ldate:
        
        print '-- Testing potential GLOBAL offset at date: ',date
        
        if tmp_gts.test_offset_significance(date , lcomponent=components, mode='detrend') and \
            (tmp_gts.extract_ndates_after_date(date,10).data.shape[0]>5) and \
            (tmp_gts.extract_ndates_before_date(date,10).data.shape[0]>5):
            
            if verbose:
                print '-- adding date ',date , ' as a definitive offset date'
            lfinal_offset_dates.append(date)
            tmp_gts = tmp_gts.add_offsets_dates([date]).detrend()
            tmp_gts.offsets_dates = []


    ###########################################################
    # STEP 2 FIND SMALL OFFSETS
    ###########################################################
     
    # days used to estimate local offset
    ndays = 7


    # new threshold
    threshold=0.9
    
            
    if verbose:
        print '-- STEP #2 SMALL OFFSET SEARCH USING T-STATISTICS THRESHOLD=',threshold

    # BEGIN LOOP

    OK=True
    iteration = 0
    while OK:

        iteration = iteration + 1
        print '-- Current wrms is      = ', tmp_gts.detrend_seasonal().wrms() * 1.E3  
        print '-- Current offset list  = ', loffsets  
        print '-- Running find_offsets_t_scan with threshold = ', threshold , ' iteration #', iteration 

        tmp_gts = tmp_gts.find_offsets_t_scan(threshold=threshold,window=50,lcomponent=components,verbose=verbose,debug=False)

        # Eliminate the offsets already tested and that were not significant
        
        if lnot_signicant_offset != []:
            tmp_gts.offsets_dates = __remove_already_tested_dates(tmp_gts.offsets_dates,lnot_signicant_offset)

        # Eliminate already found offsets
        print '-- Cheking whether these suspected offsets are new'
        
        if loffsets != []:
            tmp_gts.offsets_dates = __remove_already_tested_dates(tmp_gts.offsets_dates,loffsets)

        print '-- There are ',len(tmp_gts.offsets_dates),' potential offsets to be tested'

        # IF NO OFFSET SUSPECTED GO TO NEXT STEP             
        if tmp_gts.offsets_dates == []:
            if verbose:
                print '-- No more offset suspected for threshold = ',threshold,'. Going to next step.'
            OK=False
            break
        
        # TEST SUSPECTED OFFSETS
        if verbose:
            print '-- found ',len( tmp_gts.offsets_dates ),' potential offsets'
            
        print '-- Testing significance of potential offsets'
        
        ldate = tmp_gts.offsets_dates

        tmp_gts.offsets_dates = []
        
        for date in ldate:
            
            print '-- Testing potential offset at date: ',date,__fmt_date(date)
            
            SIGNIFICANT=tmp_gts.test_offset_significance(date, lcomponent=components,mode='local',conf_level=95)
             
            if SIGNIFICANT:
                print '-- offset at date ',date, __fmt_date(date),' is significant'
                loffsets.append( date )
                offset_values = tmp_gts.local_offset_robust(date , ndays ,verbose=verbose)
                print '-- correcting offset at date ',date, __fmt_date(date), ' with values ', offset_values[1:4]*1.E3, 'mm'
                tmp_gts = tmp_gts.apply_offsets(offset_values.reshape(-1,7),opposite=True)
                
            else:
                print '-- Offset at date ',date,__fmt_date(date),' not significant'
                lnot_signicant_offset.append(date)
        
        # test whether the maximum iteration has been reached
        if iteration >= itermax:
            print '-- Maximum number of iteration reached: ', itermax
            OK=False

    # END LOOP
        
        
    print '-- Kept offset_dates ',loffsets


    
    if conf_level_subtle >0:
    
        ###########################################################
        # STEP 3 FIND SUBTILE OFFSETS
        ###########################################################

        # new threshold
        threshold=0.8

        if verbose:
            print '-- STEP #3 SUBTLE OFFSET SEARCH USING T-STATISTICS THRESHOLD=',threshold


        tmp_gts = self.add_offsets_dates(loffsets).detrend_seasonal().find_outliers_vondrak(threshold=5., fc=2. , verbose=True, component=components).remove_outliers()
        tmp_gts.offsets_dates = []
 
        OK=True
        iteration = 0

        # BEGIN LOOP
        
        while OK:

            iteration = iteration + 1
            print '-- Current wrms is      = ', tmp_gts.detrend_seasonal().wrms() * 1.E3  
            print '-- Current offset list  = ', loffsets  
            print '-- Running find_offsets_t_scan with threshold = ', threshold , ' iteration #', iteration 

            tmp_gts = tmp_gts.find_offsets_t_scan(threshold=threshold,window=50,lcomponent=components,verbose=verbose,debug=False)

            # Eliminate the offsets already tested and that were not significant
            
            if lnot_signicant_offset != []:
                tmp_gts.offsets_dates = __remove_already_tested_dates(tmp_gts.offsets_dates,lnot_signicant_offset)

            # Eliminate already found offsets
            print '-- Cheking whether these suspected offsets are new'
            
            if loffsets != []:
                tmp_gts.offsets_dates = __remove_already_tested_dates(tmp_gts.offsets_dates,loffsets)

            print '-- There are ',len(tmp_gts.offsets_dates),' potential offsets to be tested'

            # IF NO OFFSET SUSPECTED GO TO NEXT STEP             
            if tmp_gts.offsets_dates == []:
                if verbose:
                    print '-- No more offset suspected for threshold = ',threshold,'. Going to next step.'
                OK=False
                break
            
            # TEST SUSPECTED OFFSETS
            if verbose:
                print '-- found ',len( tmp_gts.offsets_dates ),' potential offsets'
                
            print '-- Testing significance of potential offsets'
            
            ldate = tmp_gts.offsets_dates

            tmp_gts.offsets_dates = []
            
            for date in ldate:
                
                print '-- Testing potential offset at date: ',date,__fmt_date(date)
                
                SIGNIFICANT=tmp_gts.test_offset_significance(date, lcomponent=components,mode='local',conf_level=conf_level_subtle)
                 
                if SIGNIFICANT:
                    print '-- offset at date ',date, __fmt_date(date),' is significant'
                    loffsets.append( date )
                    print '-- correcting offset at date ',date, __fmt_date(date)
                    tmp_gts = self.add_offsets_dates(loffsets).detrend_seasonal().find_outliers_vondrak(threshold=5., fc=2. , verbose=True, component=components).remove_outliers().detrend_seasonal()
                    tmp_gts.offsets_dates = []
                    
                else:
                    print '-- Offset at date ',date,__fmt_date(date),' not significant'
                    lnot_signicant_offset.append(date)

            # test whether the maximum iteration has been reached
            if iteration >= itermax:
                print '-- Maximum number of iteration reached: ', itermax
                OK=False
        
    # END LOOP
        
        
    print '-- Kept offset_dates ',loffsets

    
    
    # Find Outliers
    tmp_gts.offsets_dates = loffsets
    tmp_gts = tmp_gts.detrend_seasonal().find_outliers_vondrak(threshold=5., fc=2. , verbose=True, component=components)
    
    # Final Global test on offsets

    loutliers = tmp_gts.outliers
    ldate = tmp_gts.offsets_dates

    lfinal_offset_dates = []
    
    tmp_gts = self.copy()
    
    tmp_gts.offsets_dates = []
    tmp_gts.outliers = loutliers
    tmp_gts = tmp_gts.remove_outliers() 

    for date in ldate:
        
        print '-- Testing potential GLOBAL offset at date: ',date
        
        if tmp_gts.test_offset_significance(date , lcomponent=components, mode='detrend') and \
            (tmp_gts.extract_ndates_after_date(date,10).data.shape[0]>5) and \
            (tmp_gts.extract_ndates_before_date(date,10).data.shape[0]>5):
            
            if verbose:
                print '-- adding date ',date , __fmt_date(date), ' as a definitive offset date'
            lfinal_offset_dates.append(date)
            tmp_gts = tmp_gts.add_offsets_dates([date]).detrend()
            tmp_gts.offsets_dates = []





# RETURN
#####################################
    
    print '-- FINAL RESULTS '
    print '-- NUMBER OF OFFSETS FOUND: ',len(lfinal_offset_dates)
    for date in lfinal_offset_dates:
        print("%10.5lf %s" % (date,__fmt_date(date)))
    print '-- OUTLIERS FLAGGED: ',len(loutliers)
    
    if in_place:
#            self.outliers=new_Gts.outliers
        self.offsets_dates = lfinal_offset_dates
        self.outliers = loutliers
        return(self)
#            del new_Gts
    else:
        new_Gts = self.copy()
        new_Gts.offsets_dates = lfinal_offset_dates
        new_Gts.outliers = loutliers
        return(new_Gts)
    
