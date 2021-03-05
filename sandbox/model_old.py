
import numpy as np
from pyacs.gts.Gts import Gts
import inspect

###################################################################
def add_vel_sigma(self,in_place=False, b_fn=4, verbose=True):
###################################################################
    """
    calculates realistic sigma on velocity components assuming white &  
    flicker using eq (19) & (23) from Williams (J. of Geodesy, 2003)
    b_fn is the value for flicker noise, taken as 4 mm/yr^1/4
    model can be detrend, detrend_annual, detrend_seasonal 
    if in_place = True then replace the current time series
    """
    
    # white noise
    # white noise is estimated from the time series
    
    if not isinstance(self.velocity,np.ndarray):
        print("!!!ERROR: Can't estimate velocity sigma before velocity because a residual time series is required to estimate noise components")
        print("!!!ERROR: use detrend, detrend_annual or detrend_seasonal first.")
        return(None)
    
    n=self.data.shape[0]
    if n>5:
        (a_n,a_e,a_u)=np.std(np.diff(self.data,n=1,axis=0),axis=0)[1:4]
    else:
        (a_n,a_e,a_u)=(2.0,2.0,5.0)
        
    def sigma2_wn(ti,a):
        ti=ti-ti[0]
        n=ti.shape[0]
        sigma2_wn=n*a**2 / (n*np.sum(ti**2) - (np.sum(ti))**2)
        return(sigma2_wn)
    # flicker noise
    def sigma2_fn(ti,b_fn):

        n=ti.shape[0]
        delta_t=np.mean(np.diff(ti,n=1))
        
        sigma2_fn= 9*b_fn**2 / (16 * delta_t**2 * (n**2 -1))
#        print 'duration ',delta_t**2 * n**2
        return(sigma2_fn)
    
    
    ti=self.data[:,0]

    if verbose:
        print('fn: ',sigma2_fn(ti,b_fn))
    
    new_ts=self.copy()
    min_rms=np.min(np.array([a_n,a_e,a_e]))

    sigma_vn=np.sqrt(sigma2_wn(ti,a_n)+sigma2_fn(ti,b_fn))*(a_n/min_rms)
    sigma_ve=np.sqrt(sigma2_wn(ti,a_e)+sigma2_fn(ti,b_fn))*(a_e/min_rms)
    sigma_vu=np.sqrt(sigma2_wn(ti,a_u)+sigma2_fn(ti,b_fn))*(a_u/min_rms)
    # to meters
    new_ts.velocity[3]=sigma_vn/1000.0
    new_ts.velocity[4]=sigma_ve/1000.0
    new_ts.velocity[5]=sigma_vu/1000.0

    if in_place:
        self.data=new_ts.data.copy()
    return(new_ts)


###################################################################
def detrend(self, method='L2' , in_place=False, periods=[], exclude_periods=[]):
###################################################################
    """
    detrends a time series and save velocity estimates in velocity attribute

    :param periods         : periods used to estimate the velocity
    :param exclude_periods : periods to be excluded for the velocity estimate
    :param in_place        : if True then replace the current time series
    :return                : the detrended time series
    :note                  : outliers from Gts.outliers are omitted in the estimation and
    offsets given Gts.offsets_dates are estimated simultaneously

    """
    ###########################################################################
    # check data is not None
    from pyacs.gts.lib.errors import GtsInputDataNone
    
    try:
        if self.data is None:
            # raise exception
            raise GtsInputDataNone(inspect.stack()[0][3],__name__,self)
    except GtsInputDataNone as error:
        # print PYACS WARNING
        print( error )
        return( self )
    ###########################################################################


    import copy
    outliers=copy.deepcopy(self.outliers)
    tmp_ts=self.remove_outliers()
    
    if periods != []:
        tmp_ts=tmp_ts.extract_periods(periods)
    if exclude_periods != []:
        tmp_ts=tmp_ts.exclude_periods(periods)

    detrended=tmp_ts.make_model(option='detrend',method=method)
    
    vel    =detrended.velocity
    offsets_values=detrended.offsets_values
            
    new_gts=self.copy()
    new_gts.outliers=outliers
    new_gts.offsets_values=offsets_values
    new_gts.velocity=vel
    
    model=new_gts.mmodel()
    
    new_gts.data[:,1:4]=new_gts.data[:,1:4]-model.data[:,1:4]
        
    if in_place:
            self.data=new_gts.data
            return(self)
    else:
        return(new_gts)

###################################################################
## detrend_annual
###################################################################

def detrend_annual(self,method='L2',in_place=False,periods=None,exclude_periods=None):
    """
    estimates a trend + annual terms in a time series and removes them
    velocity and annual attribute are saved in Gts.velocity & Gts.annual

    :param periods         : periods used for estimation
    :param exclude_periods : periods to be excluded from estimation
    :param in_place        : if True then replace the current time series
    :return                : the detrended time series
    :note                  : outliers from Gts.outliers are ommitted in the estimation and
    offsets given Gts.offsets_dates are estimated simultaneously

    """
    ###########################################################################
    # check data is not None
    from pyacs.gts.lib.errors import GtsInputDataNone
    
    try:
        if self.data is None:
            # raise exception
            raise GtsInputDataNone(inspect.stack()[0][3],__name__,self)
    except GtsInputDataNone as error:
        # print PYACS WARNING
        print( error )
        return( self )
    ###########################################################################

    import copy
    outliers=copy.deepcopy(self.outliers)
    tmp_ts=self.remove_outliers()
    
    if periods:
        tmp_ts=tmp_ts.extract_periods(periods)
    if exclude_periods:
        tmp_ts=tmp_ts.exclude_periods(periods)

    detrended=tmp_ts.make_model(option='detrend_annual',method=method)
    
    vel    =detrended.velocity
    offsets_values=detrended.offsets_values
    annual=detrended.annual
            
    new_gts=self.copy()
    new_gts.outliers=outliers
    new_gts.offsets_values=offsets_values
    new_gts.velocity=vel
    new_gts.annual=annual
    
    model=new_gts.mmodel()
    
    new_gts.data[:,1:4]=new_gts.data[:,1:4]-model.data[:,1:4]
        
    if in_place:
            self.data=new_gts.data
    else:
        return(new_gts)


###################################################################
## detrend_seasonal
###################################################################

def detrend_seasonal(self,method='L2',in_place=False,periods=None,exclude_periods=None):
    """
    estimates a trend + annual + semi-annual terms in a time series and removes them
    velocity, annual and semi-annual attributes are saved in Gts.velocity, Gts.annual, Gts.semi_annual

    :param periods         : periods used for estimation
    :param exclude_periods : periods to be excluded from estimation
    :param in_place        : if True then replace the current time series
    :return                : the detrended time series
    :note                  : outliers from Gts.outliers are ommitted in the estimation and
    offsets given Gts.offsets_dates are estimated simultaneously
    """

    ###########################################################################
    # check data is not None
    from pyacs.gts.lib.errors import GtsInputDataNone
    
    try:
        if self.data is None:
            # raise exception
            raise GtsInputDataNone(inspect.stack()[0][3],__name__,self)
    except GtsInputDataNone as error:
        # print PYACS WARNING
        print( error )
        return( self )
    ###########################################################################

    import copy
    outliers=copy.deepcopy(self.outliers)
    tmp_ts=self.remove_outliers()
    
    if periods:
        tmp_ts=tmp_ts.extract_periods(periods)
    if exclude_periods:
        tmp_ts=tmp_ts.exclude_periods(periods)

    detrended=tmp_ts.make_model(option='detrend_seasonal',method=method)
    
    vel    =detrended.velocity
    offsets_values=detrended.offsets_values
    annual=detrended.annual
    semi_annual=detrended.semi_annual
            
    new_gts=self.copy()
    new_gts.outliers=outliers
    new_gts.offsets_values=offsets_values
    new_gts.velocity=vel
    new_gts.annual=annual
    new_gts.semi_annual=semi_annual
    
    model=new_gts.mmodel()
    
    new_gts.data[:,1:4]=new_gts.data[:,1:4]-model.data[:,1:4]
        
    if in_place:
            self.data=new_gts.data
    else:
        return(new_gts)


###################################################################
## remove_pole
###################################################################

def remove_pole(self,pole,pole_type='euler',in_place=False, verbose=True):
    """
    remove velocity predicted by an Euler pole or a rotation rate vector from a time series
    pole is a 1D array with 3 values
    requires self.lon & self.lat attributes to have been filled before 
    if in_place = True then replace the current time series
    """

    ###########################################################################
    # check data is not None
    from pyacs.gts.lib.errors import GtsInputDataNone
    
    try:
        if self.data is None:
            # raise exception
            raise GtsInputDataNone(inspect.stack()[0][3],__name__,self)
    except GtsInputDataNone as error:
        # print PYACS WARNING
        print( error )
        return( self )
    ###########################################################################

    if (self.lon is None) or (self.lat is None):
        print("!!! ERROR: lon & lat needs to be properly filled to use this method.")
        return()
    
    from pyacs.lib.gmtpoint import GMT_Point
    
    M=GMT_Point(code=self.code,lon=self.lon,lat=self.lat,Ve=0.,Vn=0.,SVe=0.,SVn=0.)
    #N=M.substract_pole(pole,pole_type)

    N=M.pole(W=pole,SW=None,type_euler='euler',option='predict')

    vel_neu=np.array([N.Vn,N.Ve,0.])*1E-3
    
    if verbose:print("-- Removing velocity (NEU)",vel_neu*1.E3)
    
    new_Gts=self.remove_velocity(vel_neu)

    if in_place:
        self.data=new_Gts.data.copy()
    return(new_Gts)
    

###################################################################
## change reference frame for the time series
###################################################################

def frame(self,frame=None, in_place=False, verbose=False):
    """
    Rotates a time series according to an Euler pole
    Returns a new Gts instance
    """

    ###########################################################################
    # check data is not None
    from pyacs.gts.lib.errors import GtsInputDataNone
    
    try:
        if self.data is None:
            # raise exception
            raise GtsInputDataNone(inspect.stack()[0][3],__name__,self)
    except GtsInputDataNone as error:
        # print PYACS WARNING
        print( error )
        return( self )
    ###########################################################################

    # Euler poles taken from pygvel_pole_info.py
    lEuler={}
    lEuler['soam']=[-132.21,-18.83,0.121]
    lEuler['nas']=[-97.52,6.48,0.359]
    lEuler['nazca']=[-94.4,61.0,0.57]
    lEuler['nas_wrt_soam']=[-83.40,15.21,0.287]
    lEuler['inca_wrt_soam']=[-63.76,22.47,0.092]
    lEuler['eura']=[-98.83483039304333,54.22539546556553,0.25678223107826376] # from Altamimi et al., 2012, eura_wrt_itrf2008

    euler_vector=np.array(lEuler[frame])

    new_Gts=self.remove_pole(euler_vector, verbose=verbose)
    if in_place:
        self.data=new_Gts.data.copy()
    return(new_Gts)

###################################################################
##  MODEL A TIME SERIES WITH A SPLINE
###################################################################

def spline(self,smoothing=1,degree=5, date=None):
    """
    :param smoothing: Positive smoothing factor used to choose the number of knots. Number of knots will be increased
    until the smoothing condition is satisfied:
        sum((w[i] * (y[i]-spl(x[i])))**2, axis=0) <= s
    :param degree: Degree of the smoothing spline. Must be <= 5. Default is k=3, a cubic spline.
    :param date: 1D array of interpolation dates in decimal year, or 'day' for every day. defualt None will interpolate
    at data date only.
    :return: new gts instance
    """

    ###########################################################################
    # check data is not None
    from pyacs.gts.lib.errors import GtsInputDataNone
    
    try:
        if self.data is None:
            # raise exception
            raise GtsInputDataNone(inspect.stack()[0][3],__name__,self)
    except GtsInputDataNone as error:
        # print PYACS WARNING
        print( error )
        return( self )
    ###########################################################################

    # import
    from scipy.interpolate import UnivariateSpline
    import numpy as np

    Spline_ts=self.copy()
    s_e = UnivariateSpline(self.data[:,0],self.data[:,1], s=smoothing*1E-3,k=degree)
    s_n = UnivariateSpline(self.data[:,0],self.data[:,2], s=smoothing*1E-3,k=degree)
    s_u = UnivariateSpline(self.data[:,0],self.data[:,3], s=smoothing*1E-3,k=degree)
    if date is None:
        Spline_ts.data[:,1]=s_e(self.data[:,0])
        Spline_ts.data[:,2]=s_n(self.data[:,0])
        Spline_ts.data[:,3]=s_u(self.data[:,0])
    if isinstance(date,np.ndarray):
        Spline_ts.data[:,1]=s_e( date )
        Spline_ts.data[:,2]=s_n( date )
        Spline_ts.data[:,3]=s_u( date )
    if date=='day':
        import pyacs.lib.astrotime as at
        np_date = at.mjd2decyear( np.arange( at.decyear2mjd( self.data[0,0] ) , at.decyear2mjd( self.data[-1,0] ) ) )
        Spline_ts.data = np.zeros( ( np_date.shape[0] , 10 ) )
        Spline_ts.data[:,0]= np_date
        Spline_ts.data[:,1]=s_e( np_date )
        Spline_ts.data[:,2]=s_n( np_date )
        Spline_ts.data[:,3]=s_u( np_date )

    return(Spline_ts)

###################################################################
## smooth a Gts time series
###################################################################

def smooth(self,window_len=11,window='hanning', in_place=False,verbose=False, component='NEU'):
    """
    smooth a time series
    """

    ###########################################################################
    # check data is not None
    from pyacs.gts.lib.errors import GtsInputDataNone
    
    try:
        if self.data is None:
            # raise exception
            raise GtsInputDataNone(inspect.stack()[0][3],__name__,self)
    except GtsInputDataNone as error:
        # print PYACS WARNING
        print( error )
        return( self )
    ###########################################################################


    ###################################################################
    ## smoothing routines from http://wiki.scipy.org/Cookbook/SignalSmooth
    # changes numpy to np # JMN 18/07/2014
    ###################################################################
    
    def smooth_scipy(x,window_len=11,window='hanning'):
        """smooth the data using a window with requested size.
        
        This method is based on the convolution of a scaled window with the signal.
        The signal is prepared by introducing reflected copies of the signal 
        (with the window size) in both ends so that transient parts are minimized
        in the begining and end part of the output signal.
        
        input:
            x: the input signal 
            window_len: the dimension of the smoothing window; should be an odd integer
            window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
                flat window will produce a moving average smoothing.
    
        output:
            the smoothed signal
            
        example:
    
        t=linspace(-2,2,0.1)
        x=sin(t)+randn(len(t))*0.1
        y=smooth(x)
        
        see also: 
        
        numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
        scipy.signal.lfilter
     
        TODO: the window parameter could be the window itself if an array instead of a string
        NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
        """
    
        if x.ndim != 1:
            raise ValueError("smooth only accepts 1 dimension arrays.")
    
        if x.size < window_len:
            raise ValueError("Input vector needs to be bigger than window size.")
    
    
        if window_len<3:
            return x
    
    
        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
            raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")
    
    
        s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
        #print(len(s))
        if window == 'flat': #moving average
            w=np.ones(window_len,'d')
        else:
            w=eval('np.'+window+'(window_len)')
    
        y=np.convolve(w/w.sum(),s,mode='valid')
        return y



    new_east=smooth_scipy(self.data[:,1],window_len=window_len,window=window)
    new_north=smooth_scipy(self.data[:,2],window_len=window_len,window=window)
    new_up=smooth_scipy(self.data[:,3],window_len=window_len,window=window)
    
    new_Gts=self.copy()
    
    if in_place:
        return(self)
        del new_Gts
    else:
        new_Gts.data[:,1]=new_east[window_len // 2-1:new_Gts.data[:,1].shape[0]+window_len // 2-1]
        new_Gts.data[:,2]=new_north[window_len // 2-1:new_Gts.data[:,1].shape[0]+window_len // 2-1]
        new_Gts.data[:,3]=new_up[window_len // 2-1:new_Gts.data[:,1].shape[0]+window_len // 2-1]
        return(new_Gts)



######################################################################################################################
## estimate the offsets of time serie: y = a + b.t + [c.sin(2.pi.t)+d.cos(2.pi.t)]+[e.sin(4.pi.t)+f.cos(4.pi.t)] + offsets
######################################################################################################################
def make_model(self,option='detrend',method='L2',loutlier=None, in_place=False):

    """ 
        Estimate linear model parameters using least squares
        input: data: Gts format
        option are: 'detrend'/'detrend_annual'/'detrend_seasonal'
        output: new Gts object: time series is now the residuals wrt to the model and its associated values (vel, annual, semi-annual etc)
    """

    ###########################################################################
    # check data is not None
    from pyacs.gts.lib.errors import GtsInputDataNone
    
    try:
        if self.data is None:
            # raise exception
            raise GtsInputDataNone(inspect.stack()[0][3],__name__,self)
    except GtsInputDataNone as error:
        # print PYACS WARNING
        print( error )
        return( self )
    ###########################################################################

    import pyacs.lib.glinalg as glinalg
    import pyacs.lib.robustestimators as RobEst
    
#    from gts_estimators import least_square
    
    data=np.copy(self.data)
    
    ### REMOVE OUTLIERS FOR LINEAR ESTIMATION
    if loutlier: data = np.delete(data,loutlier,axis=0)
    
    ### REMOVE OFFSET DATES OUTSIDE THE TIME SPAN OF THE TIME SERIES
    
    if self.offsets_dates is not None and len(self.offsets_dates)>0:
        # keep offsets_dates only within the time series
        sel_offsets_dates=\
        [self.offsets_dates[i] for i in range(len(self.offsets_dates)) if ( self.offsets_dates[i]>self.data[0,0] and self.offsets_dates[i]<self.data[-1,0])]
        noffset = len(sel_offsets_dates)
        
        ### offsets_values = [time_offsets offsets_NEU s_offsets_NEU]
        offsets_values = np.zeros((noffset,7))
        offsets_values[:,0] = self.offsets_dates
    else:
        noffset = 0
        offsets_values = None
    
    # REF DATES
    
    t_ref=self.data[0,0]
    t_ref_seasonal=2010.0
    
    
    # INDEX
    if option =='detrend_seasonal':
        del_index = []
        n_default_unknown = 6
    elif option =='detrend_annual':
        del_index = [4,5]
        n_default_unknown = 4
    elif option =='detrend':
        del_index = [2,3,4,5]
        n_default_unknown = 2
    else:
        print('    ERROR!!! check the option of estimation: detrend/detrend_seasonal/detrend_annual')
    
    # INIT VEL, ANNUAL, SEMI_ANNUAL, RESIDUALS ARRAYS
    ### vel = [vel_N vel_E vel_U svel_N s_vel_E svel_U]
    vel = np.zeros(6)
    ### annual = [amplitude_NEU phase_NEU]
    annual = []
    ### semi_annual = [amplitude_NEU phase_NEU]
    semi_annual = []
    
    residuals = np.zeros(data.shape)
    residuals[:,0] = data[:,0]
    
    ndate = len(data[:,0])
    
    # BUILD LINEAR SYSTEM
    
    for k in range(1,4):
        
        ## write matrix A in general case
        A = np.zeros([ndate,(6+noffset)], float)
        for i in range(ndate):
            ti = data[i,0]
            A[i,0],A[i,1],A[i,2],A[i,3],A[i,4],A[i,5] = 1., (ti-t_ref),\
                                                        np.cos(2.*np.pi*(ti-t_ref_seasonal)), np.sin(2.*np.pi*(ti-t_ref_seasonal)),\
                                                        np.cos(4.*np.pi*(ti-t_ref_seasonal)), np.sin(4.*np.pi*(ti-t_ref_seasonal))
            ## for offsets
            for j in range(noffset):
                if ti > self.offsets_dates[j]: A[i,(6+j)] = 1.

        
        ### take the design matrix
        A = np.delete(A,del_index,axis=1)
        
        # solve
        
        (X,COV,V)=glinalg.lsw_full(A,data[:,k],data[:,k+3])
                
        if method == 'L1':
            (X,V)=RobEst.Dikin(A, data[:,k],data[:,k+3], eps=1.0E-4)

        
        s_X=np.sqrt(np.diag(COV))
        
        ## calculate residuals, predicted mb_files
        residuals[:,k] = V
        residuals[:,k+3] = data[:,k+3]
        
        ## velocity
        vel[k-1] = X[1]
        vel[k+2] = s_X[1]
                   
        ## calculate the offset amplitudes
        if noffset > 0:
            offsets_values[:,k] = X[n_default_unknown:(n_default_unknown+noffset)]
            offsets_values[:,k+3]=s_X[n_default_unknown:(n_default_unknown+noffset)]
            
        if option =='detrend_annual':
            ## calculate the annual motion: amplitude & phase
            annual.append((X[2],X[3]))
                            
        if option =='detrend_seasonal':
            ## calculate the annual motion: amplitude & phase
            annual.append((X[2],X[3]))
            ## calculate the annual motion: amplitude & phase
            semi_annual.append((X[4],X[5]))


    if len(annual) > 0: annual = np.hstack(annual)
    else: annual = None
    if len(semi_annual) > 0: semi_annual = np.hstack(semi_annual)
    else: semi_annual = None

    new_Gts=self.copy()
    
    new_Gts.offsets_values = offsets_values
    new_Gts.annual = annual
    new_Gts.semi_annual = semi_annual
    new_Gts.velocity = vel
    
    new_Gts.data=residuals
    
    if in_place:
        self=new_Gts.copy()
        del new_Gts
        return(self)
    else:
        return(new_Gts)

######################################################################################################################
## Generates a model time series 
######################################################################################################################

def mmodel(self):
    """
    Generates a modeled time series from the parameters read in self
    """
    
    if isinstance(self.offsets_values,np.ndarray):
        noffset = len(self.offsets_values)
        offsets_values = np.zeros((noffset,7))
    else:
        noffset = 0
        offsets_values = None
    
    t_ref=self.data[0,0]
    t_ref_seasonal=2010.0
    
    
    ### vel = [vel_N vel_E vel_U svel_N s_vel_E svel_U]
    vel = np.zeros(6)
    ### annual = [amplitude_NEU phase_NEU]
    annual = []
    ### semi_annual = [amplitude_NEU phase_NEU]
    semi_annual = []

    data=np.copy(self.data)
    
    residuals = np.zeros(data.shape)
    residuals[:,0] = data[:,0]
    
    ndate = len(data[:,0])

    for k in range(1,4):
        ## write matrix A in general case
        A = np.zeros([ndate,(6+noffset)], float)
        for i in range(ndate):
            ti = data[i,0]
            A[i,0],A[i,1],A[i,2],A[i,3],A[i,4],A[i,5] = 1., (ti-t_ref),\
                                                        np.cos(2.*np.pi*(ti-t_ref_seasonal)), np.sin(2.*np.pi*(ti-t_ref_seasonal)),\
                                                        np.cos(4.*np.pi*(ti-t_ref_seasonal)), np.sin(4.*np.pi*(ti-t_ref_seasonal))
            ## for offsets
            for j in range(noffset):
                if ti > self.offsets_dates[j]: A[i,(6+j)] = 1.
    
        X=np.zeros([(6+noffset)], float)
        
        v=self.velocity[k-1]
        if isinstance(self.annual,np.ndarray):
            annual_p=self.annual[2*(k-1)]
            annual_q=self.annual[2*(k-1)+1]
        else:
            annual_p=0
            annual_q=0
            
        if isinstance(self.semi_annual,np.ndarray):
            semi_annual_p=self.semi_annual[2*(k-1)]
            semi_annual_q=self.semi_annual[2*(k-1)+1]
        else:
            semi_annual_p=0
            semi_annual_q=0

        X[0]=0
        X[1]=v
        X[2]=annual_p
        X[3]=annual_q
        X[4]=semi_annual_p
        X[5]=semi_annual_q

        if isinstance(self.offsets_values,np.ndarray):
            for i in range(noffset):
                #print "self.offsets_values[i,k-1] ",self.offsets_values[i,k-1]
                X[6+i]=self.offsets_values[i,k]

        else:
            for i in range(noffset):
                X[6+i]=0.0

        data[:,k]=np.dot(A,X)
    new_Gts=self.copy(data=data)
    
    return(new_Gts)

###############################################################################
def detrend_median(self , delta_day = None, in_place = False, periods=[], exclude_periods=[], verbose = False, auto=False):
###############################################################################
    """
    Calculates a velocity using the median of pair of displacements exactly separated by one year, inspired from MIDAS
    If the time series has less than a year of data, then the time series is kept untouched.
    :param delta_day: if None, it is one year, if 0 then it is the relax mode for campaign data,  dany integer is the time delta (in days) used to compute velocity.
    :param in_place: boolean, if True, in_place, if False, returns a new Gts instance (default)
    :param periods: periods (list of lists) to be included for trend calculation
    :param exclude_periods: periods (list of lists) to be excluded for trend calculation
    :param verbose: verbose mode
    :param auto: if True, then start will delta_day=None, if it fails or found less than 100 pairs then use delta_day=0, 
    if fails then use regular detrend
    :note: returns None if time series is shorter than 1 year
    """    

    ###########################################################################
    # check data is not None
    from pyacs.gts.lib.errors import GtsInputDataNone
    
    try:
        if self.data is None:
            # raise exception
            raise GtsInputDataNone(inspect.stack()[0][3],__name__,self)
    except GtsInputDataNone as error:
        # print PYACS WARNING
        print( error )
        return( self )
    ###########################################################################

    ###########################################################################
    # 1-yr threshold
    ###########################################################################

    if ( self.data[-1,0] - self.data[0,0] ) < 1.:
        if verbose:
            print("!!! WARNING: time series shorter than 1 year for side: %s" % self.code)
        return( None )


    ###########################################################################
    # run the appropriate estimator using autoreference if auto is True
    ###########################################################################
    
    if auto:
        
        tts = self.detrend_median( delta_day = None, in_place = in_place, \
                                   periods=periods, exclude_periods=exclude_periods, \
                                   verbose = verbose, auto=False)

        if tts is None:
            tts = self.detrend_median( delta_day = 0, in_place = in_place, \
                                       periods=periods, exclude_periods=exclude_periods, \
                                       verbose = verbose, auto=False)
            
        if tts is None:
            tts = self.detrend( in_place = in_place, \
                                   periods=periods, exclude_periods=exclude_periods, \
                                   verbose = verbose )

        return( tts )

    ###########################################################################
    # start main detrend_median
    ###########################################################################

    import pyacs.lib.astrotime

    tmp_ts=self.remove_outliers()

    if periods != []:
        tmp_ts=tmp_ts.extract_periods(periods)
    if exclude_periods != []:
        tmp_ts=tmp_ts.exclude_periods(periods)

    # minimum number of pair velocity estimates to consider the result
    
    min_vel_estimate = 100

        
        


    
    ###########################################################################
    # MIDAS-like method, using 1-year pair data for velocity estimate
    ###########################################################################
    
    if delta_day is None:
    
        duration_in_year = tmp_ts.data[-1,0] - tmp_ts.data[0,0]
        
        if duration_in_year >= 1.: 
        
            H_year={}
            
            # creates H_year arrays
            
            for year in np.arange( int(tmp_ts.data[0,0]) , int(tmp_ts.data[-1,0])+1 ):
                
                H_year[year] = np.zeros( (365, 3) )
                H_year[year][:,:] = np.nan
                
            # deal with dates
            
            np_mjd = pyacs.lib.astrotime.decyear2mjd(tmp_ts.data[:,0])
            (np_doy,np_year,_np_ut)=pyacs.lib.astrotime.mjd2dayno(np_mjd)
            
            np_doy = np_doy.astype(int)
            np_year = np_year.astype(int)
            
            # fill H_year arrays
            
            for i in np.arange( np_doy.shape[0] ):
                if np_doy[i] != 366:
                    H_year[ np_year[i] ][ np_doy[i]-1 , : ] = tmp_ts.data[ i , 1:4 ]
        
            # stack the velocities
            
            np_vel = np.zeros((365 * (len(H_year)-1), 3)) 
            np_vel[:,:] = np.nan
            
            i=0
            for year in sorted(H_year.keys())[0:-1]:
                np_vel[i*365:(i+1)*365] = H_year[year+1] - H_year[year]
                i=i+1
            
            # test whether there are at least 3 non-Nan numbers
            
            if np.count_nonzero(~np.isnan(np_vel)) < min_vel_estimate:
                return( None )
            
            # calculates the median velocity
            
            med_vel = np.nanmedian(np_vel, axis=0)
            
            [med_vel_north, med_vel_east, med_vel_up ] = med_vel

            # calculate uncertainty and refined value of vel using step 2 from Blewitt et al. (2016) p. 2057
            # remove Nan
            
            np_vel_north = np_vel[:,0][np.logical_not(np.isnan(np_vel[:,0]))]
            np_vel_east  = np_vel[:,1][np.logical_not(np.isnan(np_vel[:,1]))]
            np_vel_up    = np_vel[:,2][np.logical_not(np.isnan(np_vel[:,2]))]
    
            # calculates sigma 
            
            fabs_vn = np.fabs(np_vel_north - med_vel_north)
            fabs_ve = np.fabs(np_vel_east  - med_vel_east )
            fabs_vu = np.fabs(np_vel_up    - med_vel_up   )
            
            sigma_vn = 1.4826 * np.median( fabs_vn )
            sigma_ve = 1.4826 * np.median( fabs_ve )
            sigma_vu = 1.4826 * np.median( fabs_vu )
            
            # removes all vel > 1.4826 * med_vel
            
            np_vel_north_cln = np_vel_north[ np.where( fabs_vn <= 2 * sigma_vn ) ]
            np_vel_east_cln  = np_vel_east[  np.where( fabs_ve <= 2 * sigma_ve ) ]
            np_vel_up_cln    = np_vel_up[    np.where( fabs_vu <= 2 * sigma_vu ) ]
    
            # get the new vel estimates
            
            med_vn = np.median( np_vel_north_cln )
            med_ve = np.median( np_vel_east_cln )
            med_vu = np.median( np_vel_up_cln )
    
            # get the new sigma
            
            med_sigma_vn = 1.4826 * np.median( np.fabs( np_vel_north_cln - med_vn ) ) 
            med_sigma_ve = 1.4826 * np.median( np.fabs( np_vel_east_cln - med_ve ) ) 
            med_sigma_vu = 1.4826 * np.median( np.fabs( np_vel_up_cln - med_vu ) ) 
    
            # empirical factor for correlated noise
            ef =3.
    
            sigma_vn = 1.2533 * med_sigma_vn / np.sqrt( fabs_vn.shape[0] ) * ef
            sigma_ve = 1.2533 * med_sigma_ve / np.sqrt( fabs_ve.shape[0] ) * ef
            sigma_vu = 1.2533 * med_sigma_vu / np.sqrt( fabs_vu.shape[0] ) * ef
    
            # final med_vel
            
            med_vel = np.array( [med_vn, med_ve, med_vu , sigma_vn, sigma_ve, sigma_vu ] )
    
        else:
            # time series has less than a year of data
            # velocity is set to 0.
            med_vel = np.array([0.,0.,0.])
            
    # end case 1-yr
    
    
    elif ( isinstance(delta_day, int) and delta_day == 0):

    ###########################################################################
    # case campaign, close to relaxed method in MIDAS
    ###########################################################################
        
        H_year={}
        lvel = [] 
        
        # creates H_year arrays
        
        for year in np.arange( int(tmp_ts.data[0,0]) , int(tmp_ts.data[-1,0])+1 ):
            
            lindex=np.where( tmp_ts.data.astype(int)[:,0]==year )
            
#            tt = tmp_ts.extract_periods([float(year),float(year)+0.999999])
#            if tt.data is not None:
#            H_year[year] = tt.data[:,:4]
            H_year[year] = tmp_ts.data[lindex[0],:4]
            
        # removes years with no available data
        # H_year = {key:val for key, val in H_year.items() if val is not None}
        # calculates velocity values for all data pairs separated by more than one year 
#        while H_year != {}:
            
        for year_start in sorted(list(H_year.keys())):
#            np_start = H_year[year_start]

            for year_end in sorted(list(H_year.keys())):
                
                ok = False
                
                for i in np.arange( H_year[ year_end ].shape[0] ):
                    for j in np.arange( H_year[ year_start ].shape[0] ):
                        
                        # check wether the pair date includes a given offset date
                        
                        include_offset_date = False 
                        
                        for odate in self.offsets_dates:
                            
                            if ( odate >=  H_year[ year_start ][j,0] ) and ( odate <= H_year[ year_end ][i,0] ):
                                include_offset_date = True

                        if not include_offset_date:
                            
                            # displacement
                            v =  (H_year[ year_end ][i] - H_year[ year_start ][j])[1:]
                            # delta_year
                            delta_year = H_year[ year_end ][i,0] - H_year[ year_start ][j,0]
                            # append to lvel
                            if delta_year > 0.90:
                                lvel.append( v/delta_year )
                                ok = True
                                if verbose:
                                    print("-- adding %.5lf - %.5lf " % (H_year[ year_start ][j,0] , H_year[ year_end ][i,0] ) )
                                    print( i,j,v/delta_year )
                                    print(H_year[ year_end ][i] , H_year[ year_start ][j])
                # if some pairs have already been found, stop here to mitigate the effect of offsets
                #if ok:
                #    break
                        

            del H_year[year_start]

                            
        # calculates the median velocity

        np_vel = np.array( lvel )
        
        np_vel_north = np_vel[:,0]
        np_vel_east  = np_vel[:,1]
        np_vel_up    = np_vel[:,2]
        
        med_vel = np.median( np_vel, axis=0)

        [med_vel_north, med_vel_east, med_vel_up ] = med_vel        
        
        # calculates sigma 
        
        fabs_vn = np.fabs(np_vel_north - med_vel_north)
        fabs_ve = np.fabs(np_vel_east  - med_vel_east )
        fabs_vu = np.fabs(np_vel_up    - med_vel_up   )
        
        sigma_vn = 1.4826 * np.median( fabs_vn )
        sigma_ve = 1.4826 * np.median( fabs_ve )
        sigma_vu = 1.4826 * np.median( fabs_vu )
        
        # removes all vel > 1.4826 * med_vel
        
        np_vel_north_cln = np_vel_north[ np.where( fabs_vn <= 2 * sigma_vn ) ]
        np_vel_east_cln  = np_vel_east[  np.where( fabs_ve <= 2 * sigma_ve ) ]
        np_vel_up_cln    = np_vel_up[    np.where( fabs_vu <= 2 * sigma_vu ) ]

        # get the new vel estimates
        
        med_vn = np.median( np_vel_north_cln )
        med_ve = np.median( np_vel_east_cln )
        med_vu = np.median( np_vel_up_cln )

        # get the new sigma
        
        med_sigma_vn = 1.4826 * np.median( np.fabs( np_vel_north_cln - med_vn ) ) 
        med_sigma_ve = 1.4826 * np.median( np.fabs( np_vel_east_cln - med_ve ) ) 
        med_sigma_vu = 1.4826 * np.median( np.fabs( np_vel_up_cln - med_vu ) ) 

        # empirical factor for correlated noise
        ef =3.

        sigma_vn = 1.2533 * med_sigma_vn / np.sqrt( fabs_vn.shape[0] ) * ef
        sigma_ve = 1.2533 * med_sigma_ve / np.sqrt( fabs_ve.shape[0] ) * ef
        sigma_vu = 1.2533 * med_sigma_vu / np.sqrt( fabs_vu.shape[0] ) * ef

        # final med_vel
        
        med_vel = np.array( [med_vn, med_ve, med_vu , sigma_vn, sigma_ve, sigma_vu ] )
    
    else:

    ###########################################################################
    # case delta_day with integer value
    ###########################################################################

        def delta(n):
            """
            Simple delta function for integer
            """
            if not n == 0:
                return(1)
            else:
                return(0)
        
        
        # first form the day vector
        
        np_mjd_data = list(map( int,pyacs.lib.astrotime.decyear2mjd(tmp_ts.data[:,0]) ))
        # so the index are
        i_np_mjd_data = np_mjd_data - np.min(np_mjd_data)

        # the void array filled with np.nan        
        
        void = np.zeros( ( np.max(np_mjd_data)-np.min(np_mjd_data) + 1 , 3 ) ) * np.nan
        
        # I fill void with the time series at existing dates
        
        void[ i_np_mjd_data ] = tmp_ts.data[:,1:4]
        
        # now reshape the vector for easy pair differentiation
        
        # if TS = void[ : , (0,1,2) ] easy diff would be achieved by working on a array W = TS.reshape(-1,delta_day)  
        # however, this requires to complete the number of lines of TS with np.nan to be able to reshape it
        
        CTS = np.zeros( ( delta_day - np.mod( void.shape[0] , delta_day)) * delta( np.mod( void.shape[0] , delta_day ) ) ) *np.nan
        
        # init med_vel
        
        med_vel = np.array([ 0. , 0. , 0. ])
        
        to_year = 365.25 / float( delta_day )
        
        for i in [0,1,2]:
            med_vel[i] = np.nanmedian( np.diff(  np.append( void[:,i] , CTS ).reshape(-1, delta_day ).T , axis=1 ) ) * to_year

    ###########################################################################
    # return detrended time series
    ###########################################################################
    
    new_gts = self.remove_velocity( med_vel )
    new_gts.outliers= self.outliers
    new_gts.offsets_values=self. offsets_values
    new_gts.velocity = med_vel
        
#    if in_place:
#            self.velocity = new_gts.velocity
#            return( self )
#    else:
#        return( new_gts )
    
    return( self.remove_velocity( med_vel, in_place = in_place) )


###############################################################################
def detrend_seasonal_median(self , wl=11, in_place=False, verbose=False):
###############################################################################
    """
    Calculates a velocity using the median of pair of displacements exactly separated by one year, inspired from MIDAS and then removes repeating yearly signal
    If the time series has less than three years of data, then the time series is kept untouched.
    
    """    

    ###########################################################################
    # check data is not None
    from pyacs.gts.lib.errors import GtsInputDataNone
    
    try:
        if self.data is None:
            # raise exception
            raise GtsInputDataNone(inspect.stack()[0][3],__name__,self)
    except GtsInputDataNone as error:
        # print PYACS WARNING
        print( error )
        return( self )
    ###########################################################################

    import pyacs.lib.astrotime

    duration_in_year = self.data[-1,0] - self.data[0,0]
    
    if duration_in_year >= 3.: 
    
        H_year={}
        
        # creates H_year arrays
        
        for year in np.arange( int(self.data[0,0]) , int(self.data[-1,0])+1 ):
            
            H_year[year] = np.zeros((365, 3))
            H_year[year][:,:] =  np.nan
            
        # deal with dates
        
        np_mjd = pyacs.lib.astrotime.decyear2mjd(self.data[:,0])
        (np_doy,np_year,_np_ut)=pyacs.lib.astrotime.mjd2dayno(np_mjd)
        
        np_doy = np_doy.astype(int)
        np_year = np_year.astype(int)
        
        # fill H_year arrays
        
        for i in np.arange( np_doy.shape[0] ):
            if np_doy[i] != 366:
    
                H_year[ np_year[i] ][ np_doy[i]-1 , : ] = self.data[ i , 1:4 ]
    
    
        # stack the velocities
        
        np_vel = np.zeros((365 * (len(H_year)-1), 3))
        np_vel[:,:] = np.nan
        
        i=0
        for year in sorted(H_year.keys())[0:-1]:
            
            np_vel[i*365:(i+1)*365] = H_year[year+1] - H_year[year]
            i=i+1
        # calculates the median velocity
        
        med_vel = np.nanmedian(np_vel, axis=0)

    
        # return detrended time series
    
        detrended = self.remove_velocity( med_vel )
    
        H_year={}
        
        # creates H_year arrays
        
        for year in np.arange( int(detrended.data[0,0]) , int(detrended.data[-1,0])+1 ):
            
            H_year[year] = np.zeros((365, 3))
            H_year[year][:,:] = np.nan
            
        # deal with dates
        
        np_mjd = pyacs.lib.astrotime.decyear2mjd(detrended.data[:,0])
        (np_doy,np_year, _np_ut)=pyacs.lib.astrotime.mjd2dayno(np_mjd)
        
        np_doy = np_doy.astype(int)
        np_year = np_year.astype(int)
        
        # fill H_year arrays
        
        for i in np.arange( np_doy.shape[0] ):
            if np_doy[i] != 366:
    
                H_year[ np_year[i] ][ np_doy[i]-1 , : ] = detrended.data[ i , 1:4 ]
    
        # center all H_year arrays
        
        for year in sorted( H_year.keys() ):
            H_year[ year ] = H_year[ year ] - np.nanmedian( H_year[ year ] , axis=0 )
#            plt.plot(H_year[year][:,1])
        # create the median daily signal

        A=np.array(list(H_year.values()))
        
        np_doy_median_signal = np.nanmedian(A[:,:,:],axis=0)

#        plt.plot(np_doy_median_signal[:,1],'ro')

        # run a median filter on it
        
        np_doy_median_signal_3 = np.vstack( (np_doy_median_signal,np_doy_median_signal,np_doy_median_signal) )
        
        import scipy.signal
        np_doy_median_signal[:,0] = scipy.signal.medfilt(np_doy_median_signal_3[:,0],wl)[365:2*365] 
        np_doy_median_signal[:,1] = scipy.signal.medfilt(np_doy_median_signal_3[:,1],wl)[365:2*365] 
        np_doy_median_signal[:,2] = scipy.signal.medfilt(np_doy_median_signal_3[:,2],wl)[365:2*365] 

#        plt.plot(np_doy_median_signal[:,1],'bo')

        # remove it from the detrended time series
        
        detrended_seasonal = detrended.copy()
        
            # loop on np_doy

        for i in np.arange( detrended_seasonal.data.shape[0] ):
            if np_doy[i] != 366:
                detrended_seasonal.data[ i , 1:4 ] = detrended_seasonal.data[ i , 1:4 ] - np_doy_median_signal[ np_doy[i]-1 , : ] 
            
    
        
    else:
        # time series is shorter than minimum

        detrended_seasonal = self.copy()

    
    return( detrended_seasonal )
