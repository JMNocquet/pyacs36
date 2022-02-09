###################################################################
def stat_site(self,lsite=[],lsite_file=None,verbose=False, save_dir=None):
###################################################################
    """
    basic statistics about individual time series
    
    :param lsite: list of selected sites for statistics 
    :param lsite_file: list of selected sites for statistics provided as a file
    :param verbose: verbose mode
    :param save_dir: directory where statistics files will be written 
     
    """
    # import
    import numpy as np
    import os , sys
     
    # Case lsite_file option
     
    if lsite_file:
        f=open(lsite_file,'r')
        for line in f:
            lline=line.split()
            if (len(lline[0])!=4):continue
            lsite.append(line.split()[0])
        f.close()
     
     
    # save_dir case
    
    if save_dir is not None:
    
        if os.path.isfile(save_dir):
            print('!!! ',save_dir,' file already exists')
            sys.exit()
         
        else:
            if not os.path.isdir(save_dir):
                try:
                    os.mkdir(save_dir)
                except:
                    print('!!! ERROR: Could not create ' , save_dir)
                    sys.exit()
    
        info_file=open(save_dir+'/pyacs.info','w')
    
     
    # header & print formats
     
    header   = "site     long.      lat.     v_e     v_n     v_u    sv_e    sv_n    sv_u   #obs #camp   s_year   e_year d_year     se     sn     su   wrms_e   wrms_n   wrms_u"
    sep_line = "--------------------------------------------------------------------------------------------------------------------------------------------------------------"
    fmt      = "%4s %9.4lf  %8.4lf %7.2lf %7.2lf %7.2lf %7.2lf %7.2lf %7.2lf %6d %5d %8.3lf %8.3lf  %5.2lf %6.2lf %6.2lf %6.2lf %8.2lf %8.2lf %8.2lf"
    fmt_vna  = "%4s %9.4lf  %8.4lf     N/A     N/A     N/A     N/A     N/A     N/A %6d %5d %8.3lf %8.3lf  %5.2lf %6.2lf %6.2lf %6.2lf %8.2lf %8.2lf %8.2lf"
    fmt_na   = "%4s %9.4lf  %8.4lf     N/A     N/A     N/A     N/A     N/A     N/A %6d %5d %8.3lf %8.3lf  %5.2lf    N/A    N/A    N/A      N/A      N/A      N/A"
     
    if save_dir:
        info_file.write(header+"\n")
        info_file.write(sep_line+"\n")
     
    print(header)
    print(sep_line)
    # loop on Gts
     
    for gts in self.lGts():
         
        # if no data move to next gts
        if gts.data is None:
            continue
         
        if lsite != [] and gts.code not in lsite:continue
        code=gts.code
        lon=gts.lon
        lat=gts.lat
        n_obs=gts.data.shape[0]
        start_year=np.min(gts.data[:,0])
        end_year=np.max(gts.data[:,0])
        delta_year=end_year-start_year
        n_camp=np.unique(list(map(int,list(gts.data[:,0])))).shape[0]
         
    #            if delta_year>2 and n_obs>2 and n_camp>1:
         
        # detrend
        new_gts=gts.detrend_median( auto=True )
    
        if new_gts is not None:
             
            [wrms_n,wrms_e,wrms_u]=new_gts.data[:,1:4].std(axis=0)*1000.0
            [s_n   ,   s_e,s_u   ]= np.median( np.fabs( np.diff( new_gts.data[:,1:4]*1000.0 , axis=0 ) ) , axis=0 )
            [vn,ve,vu,svn,sve,svu]=new_gts.velocity[:6]*1000.0
             
            fmt_str = (fmt % (code,lon,lat, ve,vn,vu, sve, svn, svu, n_obs,n_camp, start_year,end_year,delta_year,s_n,s_e,s_u,wrms_e,wrms_n,wrms_u))
            print(fmt_str)
             
            if save_dir is not None:
                info_file.write("%s\n" % (fmt_str))
                 
        else:
            if gts.data.shape[0] >= 2:
                # we can still get wrms and daily scatter
                [wrms_n,wrms_e,wrms_u]= gts.data[:,1:4].std(axis=0)*1000.0
                [s_n   ,   s_e,s_u   ]= np.median( np.fabs( np.diff( gts.data[:,1:4]*1000.0 , axis=0 ) ) , axis=0 )
    
                fmt_str = (fmt_vna % (code,lon,lat, n_obs,n_camp, start_year,end_year,delta_year,s_n,s_e,s_u,wrms_e,wrms_n,wrms_u))
                print(fmt_str)
    
                if save_dir is not None:
                    info_file.write("%s\n" % (fmt_str))
             
            else:
                fmt_str = (fmt_na % (code,lon,lat, n_obs,n_camp, start_year,end_year,delta_year))
                print(fmt_str)
    
                if save_dir is not None:
                    info_file.write("%s\n" % (fmt_str))
                 
    
    # write results
     
    if save_dir is not None:
    
        info_file.close()
        info = np.genfromtxt( save_dir+'/pyacs.info' , skip_header=2,dtype=str )
    
    
        header_info_ds = "site    long.       lat.\n"\
                       + "------------------------"
        np.savetxt(save_dir+'/pyacs_lphn.dat',info[:,:3] , fmt = "%4s %10s %10s" , header=header_info_ds )
         
        info_ds = info[ np.where(info[:,14]!='N/A') ]
        header_info_ds = "site daily_scatter_north (mm) daily_scatter_east (mm) daily_scatter_up (mm)\n"\
                       + "---------------------------------------------------------------------------"
        np.savetxt(save_dir+'/pyacs_daily_scatter.dat',info_ds[:, (0,14,15,16) ] , fmt = "%4s %20s %20s %20s" , header=header_info_ds)
         
        info_ds = info[ np.where(info[:,17]!='N/A') ]
        header_info_ds = "site    wrms_north       wrms_east    wrms_up (mm) \n"\
                       + "---------------------------------------------------"
        np.savetxt(save_dir+'/pyacs_wrms.dat',info_ds[:, (0,17,18,19) ] , fmt = "%4s %15s %15s %15s" , header=header_info_ds)
         
        np.savetxt(save_dir+'/pyacs_void.dat',info[:, (1,2,0) ] , fmt = "%10s %10s   0.00   0.00   0.00   0.00    0.00   %s")
    
        info_ds = info[ np.where(info[:,3]!='N/A') ]
        header_info_ds = "   long.       lat.         ve         vn        sve        svn    sven    site\n"\
                       + "------------------------------------------------------------------------------------"
        np.savetxt(save_dir+'/pyacs_vel.dat',info_ds[:, (1,2,3,4,6,7,0) ] , fmt = "%10s %10s %10s %10s %10s %10s    0.00    %s", header=header_info_ds)
    
        header_info_ds = "   long.       lat.        --       v_up       ---   sigma_v_up      ---    site\n"\
                       + "---------------------------------------------------------------------------------"
        np.savetxt(save_dir+'/pyacs_vel_up.dat',info_ds[:, (1,2,5,8,0) ] , fmt = "%10s %10s      0.00 %10s      0.00 %12s     0.00    %s" , header=header_info_ds)
     
     
    lcode=[]
    for gts in self.lGts():lcode.append(gts.code)
    for site in lsite:
        if site not in lcode:print("!!! ",site," in ",lsite_file," and not in ts")
