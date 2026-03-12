def __make_pytrf_model__(gts, pytrf_ts, noise = ['wn','fn'], fixed_correlated_noise_value=[None,None]):
    """
    Create a (simple) pytrf model from pyacs gts and pytrf ts
    """
    from pytrf.ts import model
    import numpy as np
    import pyacs.lib.astrotime as at

    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG

    # check noise argument consistency
    if 'fn' in noise and fixed_correlated_noise_value[1] is not None:
        ERROR("Cannot specify fixed correlated noise spectral index value when 'fn' is in noise", exit=True)
    if 'rw' in noise and fixed_correlated_noise_value[1] is not None:
        ERROR("Cannot specify fixed correlated noise spectral index value when 'rw' is in noise", exit=True)
    if len(noise) !=2:
        ERROR(f"Pyacs wrapper for PyTRF only allows two noise processes. Directly use the PyTRF software for advanced use. noise = {noise}", exit=True)
    if 'wn' not in noise:
        ERROR(f"Pyacs wrapper for PyTRF requires white noise to be included in the noise processes. Directly use the PyTRF software for advanced use. noise = {noise}", exit=True)
    
    # standard usage
    if fixed_correlated_noise_value == [None, None]:
            m = model(pytrf_ts,t0=at.decyear2mjd(int(gts.data[0, 0])), deg=[0,1],per=[365.25,182.625], noise = noise)
    else:
            m = model(pytrf_ts,t0=at.decyear2mjd(int(gts.data[0, 0])), deg=[0,1],per=[365.25,182.625], noise = ['wn'])
            noise.remove('wn')
            cnoise = noise[0]

            # set s2, fix_s2, a, fix_a consistently: value are not fixed unless fixed_correlated_noise_value is specified
            s2 = fixed_correlated_noise_value[0]
            fix_s2 = False if s2 is None else True 
            a = fixed_correlated_noise_value[1]
            fix_a = False if a is None else True 

            if cnoise == 'fn':
                  a = -1.0
                  fix_a = True
            elif cnoise == 'rw':
                  a = -2.0
                  fix_a = True
            VERBOSE(f"Adding power law noise with fixed spectral index a = {a} and variance s2 = {s2}")
            m.add_pl(a=a, fix_a=fix_a, s2=s2, fix_s2=fix_s2)

    

    # Add jumps if any
    if gts.offsets_dates is not None:
            list_offset_dates=at.decyear2mjd(gts.offsets_dates).tolist()
            m.add_jumps(list_offset_dates,deg=[0])

    return m

def __model_results_to_gts__(m, gts):
    """
    Convert pytrf model results to pyacs gts
    """
    import pyacs.lib.astrotime as at
    import numpy as np
    import pyacs.message.verbose_message as VERBOSE
    # pyacs attributes
    velocity = np.zeros(6)
    annual = np.zeros(6)
    semi_annual = np.zeros(6)
    offsets_dates_values = []
    idx2enu = {0:'East',1:'North',2:'Up'}
    idx_enu2idx_neu = {0:1,1:0,2:2}
    spectral_index = np.zeros(6)
    # Loop over deterministic parameters

    for d in range(m.nd):
        idx_offset = -1
        for f in m[d].f:
            for p in f.par:
            
                #if p.type == 'polynomial coefficient' and f.deg == 0:
                        #print("component %5s bias: %10.3f sigma: %10.3f"% (idx2enu[d],p.x*1.E3,p.sig*1.E3))

                if p.type == 'polynomial coefficient' and f.deg == 1:
                        #print("component %5s vel : %10.3f sigma: %10.3f"% (idx2enu[d],p.x*1.E3*365.25,p.sig*1.E3*365.25))
                        velocity[idx_enu2idx_neu[d]] = p.x*365.25
                        velocity[idx_enu2idx_neu[d]+3] = p.sig*365.25

                if p.type == 'polynomial coefficient jump' and f.deg == 0:
                        #print("component %5s offset %10.5f : %10.3f sigma: %10.3f"% (idx2enu[d],at.mjd2decyear(p.t),p.x*1.E3,p.sig*1.E3))
                        idx_offset += 1
                        if idx_offset >= len(offsets_dates_values):
                            offsets_dates_values.append([0,0,0,0,0,0,0])
                        offsets_dates_values[idx_offset][0] = at.mjd2decyear(p.t)
                        offsets_dates_values[idx_offset][idx_enu2idx_neu[d]+1] = p.x
                        offsets_dates_values[idx_offset][idx_enu2idx_neu[d]+4] = p.sig


                if p.type == 'cos amplitude' and f.per == 365.25:
                        #print("component %5s cos amplitude : %10.3f sigma: %10.3f"% (idx2enu[d],p.x*1.E3,p.sig*1.E3))
                        annual[idx_enu2idx_neu[d]*2] = p.x

                if p.type == 'sin amplitude' and f.per == 365.25:
                        #print("component %5s sin amplitude : %10.3f sigma: %10.3f"% (idx2enu[d],p.x*1.E3,p.sig*1.E3))
                        annual[idx_enu2idx_neu[d]*2+1] = p.x

                if p.type == 'cos amplitude' and f.per == 182.625:
                        #print("component %5s cos amplitude : %10.3f sigma: %10.3f"% (idx2enu[d],p.x*1.E3,p.sig*1.E3))
                        semi_annual[idx_enu2idx_neu[d]*2] = p.x

                if p.type == 'sin amplitude' and f.per == 182.625:
                        #print("component %5s sin amplitude : %10.3f sigma: %10.3f"% (idx2enu[d],p.x*1.E3,p.sig*1.E3))
                        semi_annual[idx_enu2idx_neu[d]*2+1] = p.x

    # Loop over noise parameters
    for d in range(m.nd):
        for n in m[d].n:
            for p in n.par:
                if 'PL spectral index' in p.type:
                        #print("component %5s power law noise spectral index: %10.3f sigma: %10.3f"% (idx2enu[d],p.x,p.sig))
                        spectral_index[idx_enu2idx_neu[d]] = p.x
                        spectral_index[idx_enu2idx_neu[d]+3] = p.sig

    # convert cos/sin into cos phase
    def __pq2aphi__(p,q):
        return(np.sqrt(p**2+q**2), np.arctan2(p,q)/(2.*np.pi))

    a_n,phi_n=__pq2aphi__(annual[0],annual[1])
    a_e,phi_e=__pq2aphi__(annual[2],annual[3])
    a_u,phi_u=__pq2aphi__(annual[4],annual[5])

    np_annual = np.array([a_n, a_e, a_u, phi_n, phi_e, phi_u])

    a_n,phi_n=__pq2aphi__(semi_annual[0],semi_annual[1])
    a_e,phi_e=__pq2aphi__(semi_annual[2],semi_annual[3])
    a_u,phi_u=__pq2aphi__(semi_annual[4],semi_annual[5])

    np_semi_annual = np.array([a_n, a_e, a_u, phi_n, phi_e, phi_u])



    with np.printoptions(precision=2, suppress=True):
        #print("velocity (mm/yr) : ", velocity * 1.0e3)
        #print("annual (mm) : ", annual * 1.0e3)
        #print("semi-annual (mm) : ", semi_annual * 1.0e3)
        #print("offsets_dates_values (mm) : ", np.array(offsets_dates_values) * 1.0e3)
        VERBOSE("power law noise spectral index [NEU]: %s " % str(spectral_index))

    # get residual time series
    new_gts = gts.copy()
    new_gts.data[:,1] = m[1].v
    new_gts.data[:,2] = m[0].v
    new_gts.data[:,3] = m[2].v
    # populate gts attributes
    new_gts.velocity = velocity
    new_gts.annual = np_annual
    new_gts.semi_annual = np_semi_annual
    new_gts.offsets_values = np.array(offsets_dates_values)

    return new_gts


##########################################################################################

def detrend_pytrf(self, noise=['wn','fn'], method='Nelder-Mead', fixed_correlated_noise_value=[None,None], log_dir= None):
    """Detrend using PyTRF; estimate velocity and offsets with realistic sigmas.

    Wrapper around pytrf (https://github.com/prebischung/pytrf, P. Rebischung).
    Simplified trajectory: constant velocity, offsets, annual/semi-annual,
    white + power-law noise. For more complex models use pytrf directly.
    Power-law noise can be fixed (useful for campaign data).

    Parameters
    ----------
    noise : list, optional
        Noise model: 'wn' (mandatory) plus 'fn', 'rw', or 'pl'. Default is ['wn','fn'].
    method : str, optional
        Optimization: 'Nelder-Mead', 'BFGS', 'CG', 'Newton', 'Powell'. Default is 'Nelder-Mead'.
    fixed_correlated_noise_value : list, optional
        [s2, a] for power-law; None to estimate. Default is [None, None].
    log_dir : str, optional
        Directory for pytrf yaml log; None = current dir. Default is None.

    Returns
    -------
    Gts
        Residual time series with velocity (and sigmas), offsets.

    Notes
    -----
    'Nelder-Mead' is default and often more robust than 'BFGS' despite being slower.
    """
    import os
    import time

    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG


    # test that pytrf is installed
    try:          
        import pytrf
    except ImportError:
        ERROR("pytrf is not installed, please install it to use this function", exit=True)

    # convert pyacs gts to pytrf ts
    r = self.to_pytrf()

    # create model
    m = __make_pytrf_model__(self, r, noise = noise, fixed_correlated_noise_value=fixed_correlated_noise_value)
    
    # fit model

    t0 = time.perf_counter()
    m.fit(estimator='reml', method=method)
    VERBOSE(f"Execution time using {method} method: {time.perf_counter() - t0:.3f} s")

    # get residuals and results
    new_gts = __model_results_to_gts__(m, self)

    # log results
    if log_dir is not None:
        if not os.path.exists(log_dir):
            os.makedirs(log_dir)
    else:
        log_dir = os.getcwd() + '/pytrf_logs'
        if not os.path.exists(log_dir):
            os.makedirs(log_dir)
    log_file = os.path.join(log_dir, self.code + '_pytrf.yaml')

    with open(log_file, "w") as f:
        print(m, file=f)

    VERBOSE("Model results saved to %s" % log_file)
    
    return new_gts