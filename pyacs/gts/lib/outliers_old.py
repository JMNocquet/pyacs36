
import numpy as np

from pyacs.gts.Gts import get_index_from_dates

###################################################################
## Remove outliers
###################################################################

###############################################################################
def remove_outliers(self,periods=None,in_place=False):
###############################################################################
    """
    removes outliers provided in self.outliers
    return a new Gts without the outliers
    if in_place = True then self has the outliers removed as well (in _place)
    """
    if self.outliers:
        if periods==None:
            data = np.delete(self.data,self.outliers,axis=0)
        else:
            lindex=self.lindex_in_periods(periods)
            ldelete=np.intersect1d(lindex,self.outliers)
            data = np.delete(self.data,ldelete,axis=0)
                
                
    else: 
        data = np.copy(self.data)
    
    new_Gts=self.copy()
    new_Gts.outliers=[]
    new_Gts.data = data
    
    if in_place:
        self.data=new_Gts.data.copy()
        del new_Gts
        self.outliers=[]
        
        return(self)
    else:
        return(new_Gts)

###################################################################
def find_outlier_around_date(self,date,conf_level=95,n=3, lcomponent='NE',verbose=True):
###################################################################
    """
    Find an outlier around a given date 
    returns the index of the outlier, returns [] if no outlier found
    :param date       : given date
    :param conf_level : confidence level for F_ratio test of outlier significance (default 95%%)
    :param n          : number of dates either sides of date (default n=3)
    :param lcomponent : components 'N','E','U','NE','NEU' (default 'NE')
    """
    
    if verbose:
        print(("-- Searching outlier around date %10.5lf on components %s with confidence level %6.1lf and %02d samples" % (date,lcomponent,conf_level,2*n)))
    
    #self.info()
    tmp_gts=self.detrend().remove_outliers().extract_ndates_around_date(date,n)
    nn=tmp_gts.data.shape[0]
    
    score={}
    llindex={}
    
    # F_ratio test
    ###############################################
    def f_ratio(chi_square_1,p1,chi_square_2,p2,n):
    ###############################################
        """
        returns result of a F_ratio test
        """
        F=( (chi_square_1-chi_square_2)/(p2-p1) ) / (chi_square_2 / (n-p2) )
        
        from scipy.stats import f
        return(f.cdf(F,p2-p1,n-p2))
    
    # 
    H_component={1:'North',2:'East',3:'Up'}
    # find outlier
    li=[]
    if 'N' in lcomponent:li.append(1)
    if 'E' in lcomponent:li.append(2)
    if 'U' in lcomponent:li.append(3)

    for i in sorted(li):
#        if verbose:
#            if i==1:print "  => Testing component: North"
#            if i==2:print "  => Testing component: East"
#            if i==3:print "  => Testing component: Up"
    
        index=np.where(np.abs(tmp_gts.data[:,i]-np.median(tmp_gts.data[:,i])) == np.max(np.abs(tmp_gts.data[:,i]-np.median(tmp_gts.data[:,i]))))
        
        if verbose:
            print(("-- suspected outlier at date %10.4lf on component %s " % (tmp_gts.data[index,0][0],H_component[i])))
        
        tmp_gts_no_outlier=tmp_gts.copy()
        tmp_gts_no_outlier.outliers=[index]
        tmp_gts_no_outlier.remove_outliers(in_place=True)
        
        chi_square_1=nn*np.std(tmp_gts.data[:,i])**2
        chi_square_2=(nn-1)*np.std(tmp_gts_no_outlier.data[:,i])**2
    
        score[i]=f_ratio(chi_square_1,1,chi_square_2,2,2*n)*100.0
        print(("-- probability of outlier (F_ratio) %5.2lf%% " % (score[i])))
        llindex[i]=index
    
    # make decision
    
    if np.max(list(score.values())) < conf_level:
        if verbose:print("-- No significant outlier found")
        return(self)
    else:
        # choose the outlier as the maximum probability
        component=li[0]
        current_score=score[component]
        del li[0]
        for i in li:
            if score[i]>current_score:
                current_score=score[i]
                component=i
        
        date=tmp_gts.data[llindex[component],0]
        # return the index in the original time series
        if verbose:
            print("=> Getting index for date ",date)
        returned_index=get_index_from_dates([date],self.data,tol=0.25)
        self.outliers+=returned_index
    
        return(self)
    
###############################################################################
def find_outliers_percentage(self,percentage=0.03,in_place=False,verbose=False, component='NEU', periods=None,excluded_periods=None):
###############################################################################
    """
    detrend a time series and
    ranks the residuals by increasing absolute value
    populate the outliers with the x % largest ones on each component
    """
    
    
    
    n_to_remove=int(percentage*self.data.shape[0])
    if n_to_remove < 1: n_to_remove=1

    if verbose:print("-- Removing the ",n_to_remove," largest outliers (each time series)")

    new_ts=self.detrend()
    lindex_north=list(np.argsort(new_ts.data[:,1]**2))[-n_to_remove:]
    lindex_east=list(np.argsort(new_ts.data[:,2]**2))[-n_to_remove:]
    lindex_up=list(np.argsort(new_ts.data[:,3]**2))[-n_to_remove:]

    lindex=[]
    if 'N' in component:lindex=lindex+lindex_north
    if 'E' in component:lindex=lindex+lindex_east
    if 'U' in component:lindex=lindex+lindex_up

    if periods and excluded_periods:
        print("!!!! periods and excluded_periods provided. Possible overlap not checked.")
        print("!!!! The program will run first on periods and then will exclude outliers in excluded_periods.")
    if periods:
        lkept_index=[]
        for index in lindex:
            for period in periods:
                start_date_period=period[0]
                end_date_period  =period[1]
                if self.data[index,0]>=start_date_period and self.data[index,0]<=end_date_period:
                    lkept_index.append(index)
                    break
        lindex=lkept_index

    if excluded_periods:
        lexcluded_index=[]
        for index in lindex:
            for period in excluded_periods:
                start_date_period=period[0]
                end_date_period  =period[1]
                if self.data[index,0]>=start_date_period and self.data[index,0]<=end_date_period:
                    lexcluded_index.append(index)
                    break
        lkept_index=[]
        for index in lindex:
            if index not in lexcluded_index:lkept_index.append()
        lindex=lkept_index
        


    new_Gts=self.copy()
    
    if in_place:
        self.outliers=list(set(lindex))
        return(self)
        del new_Gts
    else:
        new_Gts=self.copy()
        new_Gts.outliers=list(set(lindex))
        return(new_Gts)


###################################################################
## A very simple isolated outliers method
###################################################################

def find_outliers_simple(self,threshold=100,window_length=10,in_place=False,verbose=False, component='NEU', periods=None,excluded_periods=None):
    lindex_outlier=[]
    for i in range(0,self.data.shape[0]-window_length):
        window=self.data[i:i+window_length,1:4]*1000.
        #print 'window',window
        #print 'median ',np.median(window,axis=0)
        residuals=np.abs(window-np.median(window,axis=0))
        #print 'residuals',residuals
        median_of_residuals=np.median(residuals,axis=0)
        #print 'median of res ',median_of_residuals
        for index in range(0,residuals.shape[0]):
            #print residuals[index,0],median_of_residuals[0]
            
            if residuals[index,0]>threshold*median_of_residuals[0]:
                if (i+index) not in lindex_outlier:print('outlier at ',self.data[i+index,0],' N');lindex_outlier.append(i+index) 
            if residuals[index,1]>threshold*median_of_residuals[1]:
                if (i+index) not in lindex_outlier:print('outlier at ',self.data[i+index,0],' E');lindex_outlier.append(i+index) 
            if residuals[index,2]>threshold*median_of_residuals[2]:
                if (i+index) not in lindex_outlier:print('outlier at ',self.data[i+index,0],' U');lindex_outlier.append(i+index)

    if periods and excluded_periods:
        print("!!!! periods and excluded_periods provided. Possible overlap not checked.")
        print("!!!! The program will run first on periods and then will exclude outliers in excluded_periods.")
    if periods:
        lkept_index=[]
        for index in lindex_outlier:
            for period in periods:
                start_date_period=period[0]
                end_date_period  =period[1]
                if self.data[index,0]>=start_date_period and self.data[index,0]<=end_date_period:
                    lkept_index.append(index)
                    break
        lindex_outlier=lkept_index

    if excluded_periods:
        lexcluded_index=[]
        for index in lindex_outlier:
            for period in excluded_periods:
                start_date_period=period[0]
                end_date_period  =period[1]
                if self.data[index,0]>=start_date_period and self.data[index,0]<=end_date_period:
                    lexcluded_index.append(index)
                    break
        lkept_index=[]
        for index in lindex_outlier:
            if index not in lexcluded_index:lkept_index.append()
        lindex_outlier=lkept_index
        


    new_Gts=self.copy()
    
    if in_place:
        self.outliers=list(set(lindex_outlier))
        return(self)
        del new_Gts
    else:
        new_Gts=self.copy()
        new_Gts.outliers=list(set(lindex_outlier))
        return(new_Gts)



###################################################################
## find_outliers and offsets using differenciation
###################################################################
def find_outliers_and_offsets_through_differentiation(self,th=100):
    """
    find outliers and offsets using differenciation 
    """

    loutlier={}
    loffset={}
    new_Gts=self.differentiate()
    median=np.median(np.fabs(new_Gts.data),axis=0)
    lindex_max=np.argmax(new_Gts.data,axis=0)
    # test each component
    for component in (1,2,3):
        loutlier[component]=[]
        loffset[component]=[]
        index=lindex_max[component]
        #print 'th*median[component] ', th*median[component]
        if np.fabs(new_Gts.data[index,component])> th*median[component]:
            print("Suspect either an outlier or an offset at ",index, 'date ',new_Gts.data[index,0], ' on component ',component)
            # case beginning of the ts => This can only be an outlier
            if index==0:
                # then the outlier could be either the first or the second record of the time series
                if np.fabs(new_Gts.data[1,component]) > th*median[component]:
                    # the outlier is the second record
                    loutlier[component].append(1)
                else:
                    # the outlier is the first
                    loutlier[component].append(1)
                continue
            # case end of the ts => This can only be an outlier
            if index==len(new_Gts.data[:,0])-1:
                # then the outlier could be either the last or the ante-last record of the time series
                if np.fabs(new_Gts.data[len(new_Gts.data[:,0])-2,component]) > th*median[component]:
                    # the outlier is the second last record
                    loutlier[component].append(len(new_Gts.data[:,0])-2)
                else:
                    # the outlier is the first
                    loutlier[component].append(len(new_Gts.data[:,0])-1)
                continue
            # regular case
            # an isolated outlier should have two opposite AND large successive values in the differenciated time series
            # offset only has only one large value

            if np.fabs(new_Gts.data[index-1,component])< 2*median[component] and np.fabs(new_Gts.data[index+1,component])< 2*median[component]:
                # this is a jump
                loffset[component].append(new_Gts.data[index,0])
                print("-- Found offset on component ",component, " at ",new_Gts.data[index,0])
            else:
                # this is an outlier
                if np.fabs(new_Gts.data[index-1,0])>np.fabs(new_Gts.data[index+1,0]):
                    # sequence of index-1,index large differenciation
                    index=index-1
                    # else the sequence is index,index+1
                # we now test 2 criterions: similar differenciated value and
                # opposite signs
                test_large_difference=False
                test_opposite_sign=False
                if np.fabs(np.fabs(new_Gts.data[index,0])-np.fabs(new_Gts.data[index+1,0]))< 4*median[component]:test_large_difference=True
                if new_Gts.data[index,0]*new_Gts.data[index,0]:test_opposite_sign=True
                if test_large_difference and test_opposite_sign:
                    print("Found isolated outlier at : ",index)
                    loutlier[component].append(index)
                else:
                    print("Not sure what this is at date ",new_Gts.data[index,0])

    
    for outlier in loutlier[1]+loutlier[2]+loutlier[3]:        
        if self.outliers:
            if outlier not in self.outliers:self.outliers.append(outlier)
        else:
            self.outliers=[outlier]
    for offset in loffset[1]+loffset[2]+loffset[3]:
        if self.offsets:
            if offset not in self.offsets:self.offsets.append(offset)
        else:
            self.offsets=[offset]


##################################################################
## search the outliers by the time series of rms
##################################################################
def find_outliers_by_RMS_ts(self,ndays=7,th_detection=5,th_rejection=2):
    """
        Find index of outliers in a time series, populate self.outliers.
          - rms time series are first calculated over ndays
          - time windows are kept for further inspection if rms(t) > th_detection * median(rms(ts))
          - for each anomalous time windows, differentiate positions, find the max
          - test whether it is a true outlier (differentiated(t) > th_rejection * median(differentiated))
        output:
               None

    """
    data_rms = self.rms(ndays)
    
    loutlier = []
    for i in (1,2,3):

        loutlier_i = []
        out_rms = np.argwhere(data_rms[:,i] > th_detection*np.median(data_rms[:,i]))
        if len(out_rms) > 0: out_rms = np.hstack(out_rms)
        
        for j in out_rms:

            differentiate_ndays = np.hstack(np.fabs(np.diff(self.data[j:(j+ndays),i])))
            
            out_threshold = th_rejection*np.median(differentiate_ndays)
            
            index = np.argwhere(differentiate_ndays > out_threshold)
            
            if len(index) == 2:
                ### 1 outlier 
                if (index[0] == index[1]-1): loutlier_i.append(index[1] + j)
                #### 2 continuous outliers
                elif (index[0] == index[1]-2):
                    loutlier_i.append(index[0] + 1 + j)
                    loutlier_i.append(index[1] + j)

        if len(loutlier_i)!=0: 
            loutlier_i = np.hstack(loutlier_i)
            for i in np.unique(loutlier_i):
                if len(np.argwhere(loutlier_i == i)) == ndays-2: loutlier.append(i)
                elif len(np.argwhere(loutlier_i == i)) == ndays-3: loutlier.append(i)
    
    if self.outliers:
        for outlier in loutlier:
            if outlier not in self.outliers:self.outliers.append(outlier)
    else:
        self.outliers=loutlier

##################################################################
## search the outliers by the trendline
##################################################################
###############################################################################
def find_outliers_by_residuals(self,threshold=5,model='detrend_seasonal',component='NE',in_place=False):
###############################################################################
    """
    Find index of outliers by trendline/trendline_annual/trendline_seasonal (the complete model)
    Then the outliers are detected if their residuals are greater than th_rejection*standard_deviation
    
    output:
           Add the list of outlier index into self.outliers
    """


    try:
        detrended=self.make_model(option=model).data
    except:
        try:
            detrended=self.make_model(option='detrend').data
        except:
            detrended=self.data.copy()
    
#    print detrended.data
#    print detrended.data[:,1:4]
#    print np.diff(detrended.data[:,1:4],axis=0)
#    print np.abs(np.diff(detrended.data[:,1:4],axis=0),axis=0)
#    print np.median(np.abs(np.diff(detrended.data[:,1:4],axis=0),axis=0),axis=0)
    [median_north,median_east,median_up]=np.median(np.abs(np.diff(detrended[:,1:4],axis=0)),axis=0)

    lindex_north=[]
    lindex_east=[]
    lindex_up=[]
    
    if 'N' in component:
        lindex_north=np.where(np.abs(detrended[:,1]) > threshold*median_north)[0].tolist()
    if 'E' in component:
        lindex_east=np.where(np.abs(detrended[:,2]) > threshold*median_east)[0].tolist()
    if 'U' in component:
        lindex_up=np.where(np.abs(detrended[:,3]) > threshold*median_up)[0].tolist()

    loutliers=self.outliers+list(set(lindex_north+lindex_east+lindex_up))
    
    new_gts=self.copy()
    new_gts.outliers=loutliers

    if in_place:
            self.data=new_gts.data
    else:
        return(new_gts)


###############################################################################
def find_outliers_vondrak(self, threshold=10, fc=2., in_place=False,verbose=True, \
                                 periods=[[]],excluded_periods=[[]], component='NE'):
###############################################################################
    """
    Find outliers using a Vondrak filter
    """
    # init
    loutliers_dates = []
    lindex_north=[]
    lindex_east=[]
    lindex_up=[]
    
    # keep selected period
    tmp_ts=self.extract_periods(periods).exclude_periods(excluded_periods).detrend()
    
    # calculates vondrak filter
    vondrak_ts = tmp_ts.vondrak(fc,component=component,verbose=verbose)

    # calculates the residual time series

    residual_ts = tmp_ts.copy()
    
    residual_ts.data[:,1] = tmp_ts.data[:,1] - vondrak_ts.data[:,1]
    residual_ts.data[:,2] = tmp_ts.data[:,2] - vondrak_ts.data[:,2]
    residual_ts.data[:,3] = tmp_ts.data[:,3] - vondrak_ts.data[:,3]

    # get the median    
    [median_north,median_east,median_up]=np.median(np.abs( residual_ts.data[:,1:4]) , axis=0)
    
    # get the outliers
    if 'N' in component:
        lindex_north=np.where(np.abs(residual_ts.data[:,1]) > threshold*median_north)[0].tolist()
    if 'E' in component:
        lindex_east=np.where(np.abs(residual_ts.data[:,2]) > threshold*median_east)[0].tolist()
    if 'U' in component:
        lindex_up=np.where(np.abs(residual_ts.data[:,3]) > threshold*median_up)[0].tolist()

    loutliers=list(set(lindex_north+lindex_east+lindex_up))

    if verbose:print(("-- Outliers detection using Vondrak filter with fc=%.2lf : %03d new outliers detected" % (fc,len(loutliers)))) 

    # get the outliers dates
    loutliers_dates+=tmp_ts.data[loutliers,0].tolist()

    # get outliers index in original time series
    if loutliers != []:
        loutliers_index=get_index_from_dates(loutliers_dates,self.data,tol=0.25)

    
    else:
        loutliers_index=self.outliers
    
    # return
    new_Gts=self.copy()
    
    if in_place:
        self.outliers=loutliers_index
        return(self)
        del new_Gts
    else:
        new_Gts=self.copy()
        new_Gts.outliers=loutliers_index
        return(new_Gts)
    
    
    


###############################################################################
def find_outliers_sliding_window(self,\
                                 threshold=3,in_place=False,verbose=True, \
                                 periods=[[]],excluded_periods=[[]], component='NE',window_len=15,automatic=True):
###############################################################################
    """
    Find outliers using sliding windows
    """

    lindex_north=[]
    lindex_east=[]
    lindex_up=[]

    if self.data.shape[0] > window_len:
        
    
        itermax=5
    
        lindex_north=[]
        lindex_east=[]
        lindex_up=[]
        
    
    
        OK=True
        loutliers=[]
        loutliers_dates=[]
        i=0
    
        smooth=self.extract_periods(periods).exclude_periods(excluded_periods).smooth(window_len=window_len)
        new_ts=self.extract_periods(periods).exclude_periods(excluded_periods)
        residual_ts=self.extract_periods(periods).exclude_periods(excluded_periods)
        residual_ts.data[:,1:4]=new_ts.data[:,1:4]-smooth.data[:,1:4]
    
        diff_data=np.diff(self.data[:,1:4],n=1,axis=0)
        [median_north,median_east,median_up]=np.median(np.abs(diff_data),axis=0)
    
        
        
        while OK:
            
            if 'N' in component:
                lindex_north=np.where(np.abs(residual_ts.data[:,1]) > threshold*median_north)[0].tolist()
            if 'E' in component:
                lindex_east=np.where(np.abs(residual_ts.data[:,2]) > threshold*median_east)[0].tolist()
            if 'U' in component:
                lindex_up=np.where(np.abs(residual_ts.data[:,3]) > threshold*median_up)[0].tolist()
    
            loutliers=list(set(lindex_north+lindex_east+lindex_up))
    
            if verbose:print(("-- Outliers detection pass #%02d : %03d new outliers detected" % (i,len(loutliers)))) 
    
    
            #print loutliers_dates,new_ts.data[loutliers,0].tolist()
            loutliers_dates+=new_ts.data[loutliers,0].tolist()
    
            if loutliers==[]:OK=False
            
            i+=1    
            if i>itermax:OK=False
    
            smooth=self.extract_periods(periods).exclude_periods([[]]).smooth(window_len=window_len)

            
            new_ts.outliers=loutliers
            new_ts=new_ts.remove_outliers()
    
            smooth=new_ts.smooth(window_len=window_len)

    
            
            residual_ts=new_ts.copy()
            residual_ts.data[:,1:4]=new_ts.data[:,1:4]-smooth.data[:,1:4]
        
            diff_data=np.diff(self.data[:,1:4],n=1,axis=0)
            [median_north,median_east,median_up]=np.median(np.abs(diff_data),axis=0)
    
        
        if verbose:print("-- ",len(loutliers_dates)," outliers found")
        
        loutliers_index=get_index_from_dates(loutliers_dates,self.data,tol=0.25)


    else:
        loutliers_index=self.outliers

    new_Gts=self.copy()
    
    if in_place:
        self.outliers=loutliers_index
        return(self)
        del new_Gts
    else:
        new_Gts=self.copy()
        new_Gts.outliers=loutliers_index
        return(new_Gts)

