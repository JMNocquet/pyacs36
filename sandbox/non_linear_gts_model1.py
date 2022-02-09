
def nl_gts_fit(t , y , sy, model_type , offset_dates=[], eq_dates=[], H_constraints={}, bounds={}):
    
    import numpy as np
    
    def fwd_model( m ):
    
        print('tau: ' , (m[-1]) )
        
        i=0
        # trend
        y_trend = np.copy(t) * 0.0
        
        if 'trend' in model_type:
        
            y_trend = m[i] + m[i+1] * t
            i = i +2
        
        # seasonal terms
        y_annual      = np.copy(t) *.0
        y_semi_annual = np.copy(t) *.0
        
        if ('seasonal' in model_type) or ('annual' in model_type): 
            y_annual      = np.cos(2.*np.pi*(t-0.)/365.25)* m[i] + np.cos(2.*np.pi*(t-0.)/365.25)*m[i+1]
            i = i + 2 
        
        if ('seasonal' in model_type) or ('semi-annual' in model_type): 
            y_semi_annual = np.cos(4.*np.pi*(t-0.)/365.25)* m[i] + np.cos(4.*np.pi*(t-0.)/365.25)*m[i+1]
            i = i + 2 
        
        # offsets
        y_offset = np.copy(t) *.0
        
        if ('offset' in model_type):
        
            for date in offset_dates:
                
                h = t * 0. + 1.
                h[np.where( t < date )] = 0.
                y_offset = y_offset + h * m[i]
                i = i+1
        
        # psd
        
        y_psd = np.copy(t) *.0

        if ('psd' in model_type):
        
            for date in eq_dates:

                h = t * 0. + 1.
                h[np.where( t <= date )] = 0.
                y_offset = y_offset + h * m[i]
                
                h = t * 0. 
                lindex = np.where( t > date )
                
                y_psd[lindex] = y_psd[lindex] + m[i+1]*np.log( 1 + (t[lindex] - date)/np.sqrt(m[i+2]**2) )
                
                i = i + 3
        
        # sum everything
        
        y_model = y_trend + y_annual + y_semi_annual + y_offset + y_psd
        
        return y_model
#------------------------------

    def cost ( m ):

        y_model = fwd_model( m )
        # chi2_obs
        
        chi2_obs = np.sum( (y - y_model)**2 ) / y.shape[0]
        
        # chi2_reg
        
        chi2_reg = 0.
        #chi2_reg = 1000. * (m[-1] -1 )**2
        
        # chi2
        
        chi2 = chi2_obs + chi2_reg
    
        print('--chi2: ', chi2)
    
        return chi2
#---------------------------------------

    # prepare
    
    import numpy as np
    
    H_param={}
    n_param = 0
    
    if 'trend' in model_type:
        n_param = n_param + 2

    if 'seasonal' in model_type:
        n_param = n_param + 4
    
    if 'annual' in model_type:
        n_param = n_param + 2

    if 'semi-annual' in model_type:
        n_param = n_param + 2

    if 'offset' in model_type:
        n_param = n_param + len( offset_dates )

    if 'psd' in model_type:
        n_param = n_param + 3*len( eq_dates )
    
    from scipy.optimize import minimize
    
    print('-- number of parameters: ', n_param)
    
    m0= np.zeros(n_param) + 1.
    
    fwd_model( m0 )
    
    res = minimize( cost , m0, method='BFGS', tol=1e-6)
    
    print('-- minimization finished')
    print(res.message)
    print(res.x)

    print('-- plotting results')
    t_save = t
    
    from pyacs.lib import astrotime
    
    #t = astrotime.mjd2decyear( np.arange( astrotime.decyear2mjd(t[0]) , astrotime.decyear2mjd(t[-1]) ) )
    
    y_model = fwd_model(res.x)
    
    import matplotlib.pyplot as plt

    plt.plot(t_save,y)
    plt.plot(t,y_model)
    plt.show()
    
########################################################################################
# MAIN
########################################################################################

from pyacs.lib import astrotime as at
from pyacs.gts.Gts import Gts

raw_ts = Gts()
my_ts = raw_ts.read_pos(tsfile='/Users/nocquet/projets/2018/soam_proc/test_pyacs_061_snx/pos/IBEC.pos')

t= at.decyear2mjd(my_ts.data[:,0])
y=my_ts.data[:,2] * 1.E3
sy=my_ts.data[:,5]  * 1.E3

offset_dates=[at.decyear2mjd(2012.524590164)]
#offset_dates=[]
eq_dates=[at.decyear2mjd(2016.292349727)]
#eq_dates=[]
H_constraints={}
bounds={}

nl_gts_fit(t,y,sy,'trend-psd-seasonal',offset_dates = offset_dates, eq_dates=eq_dates, H_constraints=H_constraints, bounds=bounds)
