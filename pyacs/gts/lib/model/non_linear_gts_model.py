

###############################################################################
def nl_gts_fit(ty , ym , sym , model_type , offset_dates=[], eq_dates=[], H_fix={} , H_constraints={}, H_bounds={}, verbose=False):
###############################################################################
    """
    Fits a single 1D time series with a (non-linear) trajectory model.
    See documentation of gts.fit_trajectory for explanation on the options.
    """
    import numpy as np
    
###############################################################################
    def fwd_model( m ):
###############################################################################
        
        # index of parameter        
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
            y_annual      = np.cos(2.*np.pi*(t-0.)/365.25)* m[i] - np.sin(2.*np.pi*(t-0.)/365.25)*m[i+1]
            i = i + 2 
        
        if ('seasonal' in model_type) or ('semi-annual' in model_type): 
            y_semi_annual = np.cos(4.*np.pi*(t-0.)/365.25)* m[i] - np.sin(4.*np.pi*(t-0.)/365.25)*m[i+1]
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

        if ('psd_log' in model_type):
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

###############################################################################
    def cost ( m ):
###############################################################################

        y_model = fwd_model( m )
        # chi2_obs
        red_chi2_obs = np.sum( ((y - y_model)/sy)**2 ) / y.shape[0]
        
        # chi2_reg
        chi2_reg = 0.
        if H_constraints != {}:
            for k,v in H_c_prior_sigma.items():
                chi2_reg = chi2_reg + ( (m[ k ] - v[0] )/v[1] )**2
        
        # chi2
        chi2 = red_chi2_obs + chi2_reg
    
        return chi2

###############################################################################
    def rms ( m ):
###############################################################################

        y_model = fwd_model( m )

        return np.sqrt( np.sum( ((y - y_model)/1.)**2 ) / y.shape[0] )

###############################################################################
    def wrms ( m ):
###############################################################################

        y_model = fwd_model( m )

        return np.sqrt( np.sum( ((y - y_model)/1.)**2 ) / np.sum( 1./sy**2) )

###############################################################################
    def red_chi_square ( m ):
###############################################################################

        y_model = fwd_model( m )

        return np.sqrt( np.sum( ((y - y_model)/sy)**2 ) / y.shape[0] )


###############################################################################
# prepare: log parameters index 
###############################################################################
    
    H_param={}
    n_param = 0
    
    if 'trend' in model_type:
        H_param['trend_cst'] = n_param
        H_param['trend_slope'] = n_param + 1
        n_param = n_param + 2


    if 'seasonal' in model_type:
        H_param['annual_a'] = n_param
        H_param['annual_b'] = n_param + 1
        H_param['semi_annual_a'] = n_param + 2
        H_param['semi_annual_b'] = n_param + 3
        n_param = n_param + 4
    
    if 'annual' in model_type:
        H_param['annual_a'] = n_param
        H_param['annual_b'] = n_param + 1
        n_param = n_param + 2

    if 'semi-annual' in model_type:
        H_param['semi_annual_a'] = n_param 
        H_param['semi_annual_b'] = n_param + 1
        n_param = n_param + 2

    if 'offset' in model_type:
        for i in np.arange( len(offset_dates) ):
            H_param['offset'+("_%02d" % i)] = n_param + i 
        n_param = n_param + len( offset_dates )
        

    if 'psd_log' in model_type:
        for i in np.arange( len( eq_dates ) ):
            H_param['psd_log_offset'+("_%02d" % i)] = n_param + 3*i 
            H_param['psd_log_amp'+("_%02d" % i)] = n_param + 3*i + 1 
            H_param['psd_log_tau'+("_%02d" % i)] = n_param + 3*i + 2 
        n_param = n_param + 3*len( eq_dates )


    if verbose:
        print('-- number of parameters for trajectory model: ', n_param)
        print('-- index  of parameters for trajectory model: ')
        for k , v in H_param.items():
                print("-- m[%02d]: %s " % (v, k) )
        

###############################################################################
# dealing with constraints 
###############################################################################

    import scipy.optimize
    from pyacs.lib import astrotime as at
    
    # work with days rather than dec.year
    t  = at.decyear2mjd(ty)
    
    offset_dates = at.decyear2mjd(offset_dates)
    eq_dates = at.decyear2mjd(eq_dates)
    
    # work with mm rather than m
    y  = ym * 1.E3
    sy = sym * 1.E3
    
    bounds = None
    m0= np.zeros(n_param) + 1.


    # bounds either as fixed or bounded parameters
    
    if H_fix != {} or H_bounds != {}:
        
        lb = [-np.inf] * n_param
        ub = [np.inf] * n_param
        bounds = scipy.optimize.Bounds(lb, ub, keep_feasible=False)  
    
    
        # case of fixed parameters
        if H_fix != {}:
            if verbose:
                print('-- dealing with parameters hold fixed')

            for k , v in H_fix.items():
                
                lb[ H_param[ k ] ] = v
                ub[ H_param[ k ] ] = v
                m0[ H_param[ k ] ] = v
            
            bounds = scipy.optimize.Bounds(lb, ub, keep_feasible=False)  
    
            if verbose:
                print('-- bounds accounting for fixed parameters')
                for j in np.arange( len(lb) ):
                    print("-- m[%02d]: %.3lf %.3lf " % (j, lb[j], ub[j]) )


        # case of bounded parameters
        if H_bounds != {}:
            if verbose:
                print('-- dealing with bounded parameters')
            for k , v in H_bounds.items():
                
                lb[ H_param[ k ] ] = v[0]
                ub[ H_param[ k ] ] = v[1]
            
            bounds = scipy.optimize.Bounds(lb, ub, keep_feasible=False)  

            if verbose:
    
                print('-- bounds accounting for bounded parameters')
                for j in np.arange( len(lb) ):
                    print("-- m[%02d]: %.3lf %.3lf " % (j, lb[j], ub[j]) )
    
    # constraints
    
    if H_constraints != {}:

        if verbose:
            print('-- dealing with parameter constraints')
        H_c_prior_sigma = {}
        for k , v in H_constraints.items():
            H_c_prior_sigma[ H_param[ k ] ] = [ v[0] , v[1] ]
            m0[ H_param[ k ] ] = v[0]
            if verbose:
                print("-- m[%02d]: prior: %.3lf sigma: %.3lf " % ( H_param[ k ], v[0] , v[1] ) )
    
    
    
    fwd_model( m0 )
    
    res = scipy.optimize.minimize( cost , m0, method='SLSQP', tol=1e-6 , bounds=bounds)
    
    if verbose:
        print('-- minimization finished')
        print('-- scipy.optimize.minimize message: ',res.message)

    # print results

    H_res={}
    for k , v in H_param.items():
            H_res[k] = res.x[v]
            if verbose:
                print("-- %-18s m[%02d]: %+15.3lf " % (k, v,res.x[v]) )

    if verbose:
        # print rms
        print('--rms  (mm): %.2lf ' % rms( res.x ))
    
        # print rms
        print('--wrms (mm): %.2lf ' % wrms( res.x ))

        # print reduced chi2
        print('--reduced chi-square : %.2lf ' % red_chi_square( res.x ))

    # return estimated parameters , model at t, residuals at t , model at every day
    
    # parameters were already stored in H_res
    
    y_model = fwd_model(res.x)
    model =  np.vstack((at.mjd2decyear(t),y_model*1.E-3)).T
    
    residuals = np.vstack((at.mjd2decyear(t), (y-y_model)*1.E-3 )).T
    
    t_every_day       = np.arange( t[0] , t[-1] )
    t = t_every_day
    y_model_every_day = fwd_model(res.x)
    
    model_every_day = np.vstack( ( at.mjd2decyear(t_every_day),(y_model_every_day)*1.E-3) ).T
    
    return H_res , model, residuals, model_every_day

