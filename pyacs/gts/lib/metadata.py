###################################################################
def read_lon_lat(self,gmt_file, verbose=False):
###################################################################
    """
    Reads a gmt psvelo file and populates Gts.lon & Gts.lat
    
    :param gmt_file: gmt psvelo file
    :param verbose: verbose mode (boolean)

    :return: the current Gts instance
    """
    
    from pyacs.lib.vel_field import Velocity_Field

    if verbose:
        print("-- Reading ",gmt_file)

    vel=Velocity_Field()
    vel.read(file_name=gmt_file)
    
    M=vel.site(self.code)
    self.lon=M.lon
    self.lat=M.lat

    return(self)
    
###################################################################
def save_velocity(self,gmt_file, verbose=True, comment=None, up=False):
###################################################################
    """
    Appends velocity estimates (with uncertainties) to a gmt psvelo file
    
    :param gmt_file: output gmt psvelo file (will append if gmt_file already exists)
    :param verbose: verbose mode (boolean)
    :param comment: comment as a string. '# ' is pre-prended to comment if not provided 
    :param up: boolean. If True, then Ve, SVe and SVen are set to 0 and Vu and Vu are written as 4-th and 6-th fields
    
    :return: the current Gts instance
    """

    # comment
    if comment is not None:
    
        if comment[0] != '#':
            comment = '# '+comment
    else:
        comment=''

    # import
    from pyacs.lib.vel_field import Velocity_Field
    from pyacs.lib.gmtpoint import GMT_Point
    import os.path

    # init Velocity_Field() instance
    vel_field=Velocity_Field()
    
    if os.path.isfile(gmt_file):vel_field.read(gmt_file,verbose=False)
    [vn,ve,vu,svn,sve,svu]=self.velocity*1000.0
    
    if not up:
        M=GMT_Point(code=self.code+' '+comment,lon=self.lon,lat=self.lat,Ve=ve,Vn=vn,SVe=sve,SVn=svn,SVen=0.0)
    else:
        M=GMT_Point(code=self.code+' '+comment,lon=self.lon,lat=self.lat,Ve=0.0,Vn=vu,SVe=0.0,SVn=svu,SVen=0.0)
        
    vel_field.add_point(M)
    vel_field.write(gmt_file,verbose=False)
    
    return(self)

###################################################################
def save_offsets(self,ofile, verbose=True, comment='', up=False, info=False):
###################################################################
    """
    Appends offsets values to a given text file (gmt psvelo format)
    
    :param ofile: output offset file
    :param verbose: verbose mode (boolean)
    :param comment: comment as a string. '# ' is pre-prended to comment if not provided 
    :param up: boolean. If True, then Ve, SVe and SVen are set to 0 and Vu and Vu are written as 4-th and 6-th fields

    :return: the current Gts instance

    
    """

    if verbose:
        print('-- Appending offset values for site ',self.code,' into ',ofile)

    if comment != '':
        if comment[0] != '#':
            comment = '# '+comment

    
    import os.path
    import pyacs.lib.astrotime

    
    if os.path.isfile(ofile):
        f_offsets=open(ofile,'a')
    else:
        f_offsets=open(ofile,'w')

    for index in range(self.offsets_values.shape[0]):
        (date,dn,de,du,sdn,sde,sdu)=tuple(self.offsets_values[index,:]*1000.0)

        date=date/1000.0
        doy,_ut=pyacs.lib.astrotime.decyear2dayno(date)
        year=int(date)

        
        s_date=("%9.4lf = %04d %03d" % (date, year, doy))
        field_8='#'+self.code+' '+s_date+' '+comment+("%10.5lf %10.5lf" %(du,sdu))
        if up:
            f_offsets.write("%10.5lf  %10.5lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %s\n" % (self.lon,self.lat,0,du,0,sdu,0.0,self.code))
        else:
            if info:
                f_offsets.write("%10.5lf  %10.5lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %s\n" % (self.lon,self.lat,de,dn,sde,sdn,0.0,field_8))
            else:
                f_offsets.write("%10.5lf  %10.5lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %s\n" % (self.lon,self.lat,de,dn,sde,sdn,0.0,self.code))
            
    f_offsets.close()
    return(self)

###################################################################
def read_eq_rename(self,eq_rename,in_place=False,verbose=False):
###################################################################
    """
    Reads the information for the current site (code) from an eq_rename globk file.
    
    Populates loutliers and offsets_dates
    Found excluded periods in the eq_rename file are added to loutliers
    
    :param eq_rename: eq_rename (globk format) file to be read
    :param in_place: boolean. If True then the Gts instance is modified, if False the Gts instance is preserved and a new Gts instance is return
    :param verbose: verbose mode (boolean)

    :return: Gts instance
    """
    
    import pyacs.lib.astrotime
    import numpy as np

    new_Gts=self.copy()
    
    try:
        file_eq_rename=open(eq_rename,'r')
    except:
        print("-- Could not open ",eq_rename) 
        return(None)
    
    lcontent=file_eq_rename.readlines()
    
    lsite_content=[]
    
    for line in lcontent:
        if self.code in line:lsite_content.append(line)
        
    if verbose:
        print("-- Found ",len(lsite_content)," commands for site ",self.code," in ",eq_rename)
    
    if lsite_content!=[]:
        for line in lsite_content:
            lline=line.split()
            # offsets_dates
            if 'break' in line and self.code in line:
                print(line)
                print(lline[2:6])
                (year,month,mday,h,m)=list(map(int,lline[2:7]))
                decyear=pyacs.lib.astrotime.cal2decyear(mday, month, year, ut=h/24.+m/24./60.)
                new_Gts.offsets_dates.append(decyear)
                
                if verbose:
                    print(("-- Adding offset @ %4d %02d %02d %02d %02d = %10.5lf" %  (year,month,mday,h,m,decyear)))
                
            # exclude
            if (('_XPS' in line) or ('_XCL' in line)) and ('rename' in line) and self.code in line:
                (syear,smonth,smday,sh,sm,eyear,emonth,emday,eh,em)=list(map(int,lline[3:13]))
                sdecyear=pyacs.lib.astrotime.cal2decyear(smday, smonth, syear, ut=sh/24.+sm/24./60.)
                edecyear=pyacs.lib.astrotime.cal2decyear(emday, emonth, eyear, ut=eh/24.+em/24./60.)
                lindex_tuple=np.where((self.data[:,0] > sdecyear) & (self.data[:,0] < edecyear))
                lindex=list(lindex_tuple[0])
                if verbose:print(("-- Flagging outliers @ %4d %02d %02d %02d %02d %4d - %02d %02d %02d %02d = %10.5lf - %10.5lf" %  (syear,smonth,smday,sh,sm,eyear,emonth,emday,eh,em, sdecyear,edecyear)))
                new_Gts.outliers+=lindex
                
            # rename without '_XGPS' or '_XCL' indicate an offset
            if ('rename' in line) and (not (('_XPS' in line) or ('_XCL' in line))):
                (syear,smonth,smday,sh,sm,eyear,emonth,emday,eh,em)=list(map(int,lline[3:13]))
                sdecyear=pyacs.lib.astrotime.cal2decyear(smday, smonth, syear, ut=sh/24.+sm/24./60.)
                edecyear=pyacs.lib.astrotime.cal2decyear(emday, emonth, eyear, ut=eh/24.+em/24./60.)
                decyear=pyacs.lib.astrotime.cal2decyear(emday, emonth, eyear, ut=eh/24.+em/24./60.)
                new_Gts.offsets_dates.append(decyear)
                
    file_eq_rename.close()

    new_Gts.outliers=sorted(list(set(new_Gts.outliers)))
    new_Gts.offsets_dates=sorted(list(set(new_Gts.offsets_dates)))
    
    for date in new_Gts.offsets_dates:
        if date >=  new_Gts.data[-1,0]:new_Gts.offsets_dates.remove(date)

    if in_place:
        self.outliers=new_Gts.outliers
        self.offsets_dates=new_Gts.offsets_dates
        return(self)
    else:
        return(new_Gts)
    

###################################################################
def save_eq_rename(self,eq_rename,verbose=False, excluded_periods=None):
###################################################################
    """
    save results of a Gts analysis in globk format eq_rename
    
    :param eq_rename: output eq_rename file (Golbk format)
    :param verbose: verbose mode (boolean)
    :param exluded_periods: periods to be excluded

    :return: Gts instance
    
    
    """
    import pyacs.lib.astrotime
    import numpy as np

    import os,sys
    
    ###################################################################
    def __check_and_add_info__(new_info,eq_rename,verbose=False):
        """
        check whether the information is already in the eq_rename
        Info can be:
        ('XCL',date)
        ('XCL',period)
        ('break',date)
        """
        
        if verbose:print("-- Info to be tested:",new_info.rstrip())
        
        code=new_info.split()[1]
        
        
        PRINT_INFO=True
        
        # open or create the eq_rename file
        if not os.path.isfile(eq_rename):
            if verbose:print("-- Creating eq_rename file: ",eq_rename)
            lcontent=[]
        else:
            #if verbose:print "=> Opening existing eq_rename file:",eq_rename
            file_eq_rename=open(eq_rename,'r')
            lcontent=file_eq_rename.readlines()
            if verbose:print("=> Read ",len(lcontent)," commands in ",eq_rename)
            file_eq_rename.close()
        
        tolerance_for_break=3./365.25
    
        
        if new_info in lcontent:
            if verbose:print("-- ",new_info.rstrip()," already existing. Skipping.")
            PRINT_INFO=False
        else:
            for line_content in lcontent:
                lline_content=line_content.split()
                # offsets_dates
                if ('break' in line_content and code in line_content) and 'break' in new_info:
                    (year,month,mday,h,m)=list(map(int,lline_content[2:7]))
                    existing_decyear=pyacs.lib.astrotime.cal2decyear(mday, month, year, ut=h/24.+m/24./60.)
                    (year,month,mday,h,m)=list(map(int,new_info.split()[2:7]))
                    new_decyear=pyacs.lib.astrotime.cal2decyear(mday, month, year, ut=h/24.+m/24./60.)
                    if np.sqrt((existing_decyear-new_decyear)**2) < tolerance_for_break:
                        
                        if verbose:
                            print("-- Conflicting info ",new_info.rstrip(),' --with-- ',line_content.rstrip()," (not enough data between the two break dates)")
                        
                        new_break_date=(existing_decyear+new_decyear)/2.
                        year=int(new_break_date)
                        (mday,month,_ut)=pyacs.lib.astrotime.decyear2cal(new_break_date)
                        corrected_info=((" break  %s %4d %02d %02d 12 00\n")%(self.code,year,month,mday))
                        
                        if verbose:
                            print("-- Replacing",line_content.rstrip()," in ",eq_rename," --by-- ",corrected_info)
                        lcontent.remove(line_content)
                        lcontent.append(corrected_info)
                        PRINT_INFO=False

                              
                                            
                # exclude and outliers
                if ((('_XPS' in line_content or '_XCL' in line_content) and 'rename' in line_content) and code in line_content) and ('rename' in new_info):
                    try:
                        (syear,smonth,smday,sh,sm,eyear,emonth,emday,eh,em)=list(map(int,lline_content[3:13]))
                    except:
                        print("!!! Error ",line_content)
                        print((syear,smonth,smday,sh,sm,eyear,emonth,emday,eh,em))
                        sys.exit()
                        
                    existing_sdecyear=pyacs.lib.astrotime.cal2decyear(smday, smonth, syear, ut=sh/24.+sm/24./60.)
                    existing_edecyear=pyacs.lib.astrotime.cal2decyear(emday, emonth, eyear, ut=eh/24.+em/24./60.)
                    lindex_existing=list(np.where((self.data[:,0] > existing_sdecyear) & (self.data[:,0] < existing_edecyear))[0])

                    try:
                        (syear,smonth,smday,sh,sm,eyear,emonth,emday,eh,em)=list(map(int,new_info.split()[3:13]))
                    except:
                        print("Unexpected error:", sys.exc_info()[0])
                        print("!!! Error ",line_content)
                        print((syear,smonth,smday,sh,sm,eyear,emonth,emday,eh,em))
                        sys.exit()
                    new_sdecyear=pyacs.lib.astrotime.cal2decyear(smday, smonth, syear, ut=sh/24.+sm/24./60.)
                    new_edecyear=pyacs.lib.astrotime.cal2decyear(emday, emonth, eyear, ut=eh/24.+em/24./60.)
                    lindex_new=list(np.where((self.data[:,0] > new_sdecyear) & (self.data[:,0] < new_edecyear))[0])
                    
                    # check overlap
                    
                    def __intersect__(a, b):
                        """ return the intersection of two lists """
                        return list(set(a) & set(b))
                
                    def __union__(a, b):
                        """ return the union of two lists """
                        return list(set(a) | set(b))

                    if __intersect__(lindex_new,lindex_existing) != []:
                        if verbose:
                            print("-- Conflicting information ",new_info.rstrip()," --with-- ",line_content.rstrip())
                        new_period=[min(existing_sdecyear,new_sdecyear),max(existing_edecyear,new_edecyear)]
                        [sdate,edate]=new_period
                        syear=int(sdate)
                        (smday,smonth,_sut)=pyacs.lib.astrotime.decyear2cal(sdate)
                        eyear=int(edate)
                        (emday,emonth,_eut)=pyacs.lib.astrotime.decyear2cal(edate)
                        corrected_info=((" rename  %s %s_XCL %4d %02d %02d 00 00 %4d %02d %02d 23 59\n")%(self.code,self.code,syear,smonth,smday,eyear,emonth,emday))
                        if verbose:
                            print("-- Replacing ",line_content.rstrip()," --with-- ",corrected_info.rstrip())
                        lcontent.remove(line_content)
                        lcontent.append(corrected_info)
                        PRINT_INFO=False
                            
            
            if PRINT_INFO:
                if verbose:print("-- Added ",new_info," in ",eq_rename)
                lcontent.append(new_info)
                

            file_eq_rename=open(eq_rename,'w+') 
            file_eq_rename.writelines(["%s" % item  for item in sorted(set(lcontent))])
            file_eq_rename.close()
            
    ###################################################################
  
    
    # outliers
    for index_outlier in sorted(self.outliers):
        dec_year=self.data[index_outlier,0]
        year=int(dec_year)
        (mday,month,_ut)=pyacs.lib.astrotime.decyear2cal(dec_year)
        new_info=((" rename  %s %s_XCL %4d %02d %02d 00 00 %4d %02d %02d 23 59\n")%(self.code,self.code,year,month,mday,year,month,mday))
        
        __check_and_add_info__(new_info,eq_rename,verbose=verbose)
        
    
    # offsets
    for dec_year in sorted(self.offsets_dates):
        year=int(dec_year)
        (mday,month,_ut)=pyacs.lib.astrotime.decyear2cal(dec_year)
        new_info=((" break  %s %4d %02d %02d 12 00\n")%(self.code,year,month,mday))
        __check_and_add_info__(new_info,eq_rename,verbose=verbose)
        
    # excluded_periods
    if verbose:
        print("-- excluded_periods ",excluded_periods)
    if isinstance(excluded_periods,list):
        if verbose:print("=> Now excluding periods")
        for period in excluded_periods:
            [sdate,edate]=period
            syear=int(sdate)
            (smday,smonth,_sut)=pyacs.lib.astrotime.decyear2cal(sdate)
            eyear=int(edate)
            (emday,emonth,_eut)=pyacs.lib.astrotime.decyear2cal(edate)
            new_info=((" rename  %s %s_XCL %4d %02d %02d 00 00 %4d %02d %02d 23 59\n")%(self.code,self.code,syear,smonth,smday,eyear,emonth,emday))
            
            __check_and_add_info__(new_info,eq_rename,verbose=verbose)
            
            if verbose:
                print((("=> Adding exclusion cmd   %s in eq_rename file: %s")%(new_info,eq_rename)))
                
    return()

###################################################################
def make_dynamic_apr(self,apr, time_step=30., pos_tol=0.03, dates=[], gap=20.0, verbose=False):
###################################################################
    """
    Creates an apr file for GAMIT
    The created apr file has no velocity, but a series of coordinates at different time
    
    
    :param apr: apr file (Globk format)
    :param time_step: time step for writing dates (default 30 days)
    :param pos_tol: position tolerance. If exceeded, a new date will be written. (default 0.03 m)
    :param dates: a list of dates in decimal years where writing will be forced
    :param gap: gap in days. If there is no data during a duration greater than gap, then observation is forced to be included and tested against pos_tol
    :param verbose: verbose mode (boolean)

    :return: Gts instance
    
    """

    # import
    
    import numpy as np
    import os
    import pyacs.lib.astrotime

    str_comment_gen='from pyacs.gts '

    # find gaps and add observation
    
    all_dates = self.data[:,0]
    dates_delta = np.diff(all_dates)
    
    lindex = np.where( dates_delta  > (gap / 365.25) )[0] + 1

    if verbose:
        print('-- found ',lindex.shape[0],' gaps longer than ',gap,' days')
    
    if lindex.shape[0] > 0:
        data_gap = self.data[lindex,:]
        dates_gap = self.data[lindex,0]
    
    else:
        dates_gap = np.array([])
        data_gap = None

    
    # decimate gts
    
    decim_ts = self.decimate(time_step=time_step,dates=[],method='median',verbose=verbose)

    if data_gap is not None:
        decim_ts.data = np.vstack( (decim_ts.data,data_gap) )

    # neu2xyz and reorder
    
    decim_ts.neu2xyz()
    decim_ts.reorder()

    # keep only dates according to pos_tol
    
    current_pos = decim_ts.data[0,:]
    
    keep_data = decim_ts.data[0,:].reshape(1,-1)
    
    for i in np.arange(1,decim_ts.data.shape[0]):
        diff_pos = np.sqrt( np.sum( ( decim_ts.data[i,1:4] - current_pos[1:4] )**2  ) )
    
        if diff_pos > pos_tol:
            keep_data = np.vstack( (keep_data, decim_ts.data[i,:].reshape(1,-1) ) )
            current_pos = decim_ts.data[i,:]
    
            if verbose:
                date = decim_ts.data[i,0]
                (mday,month,_ut) = pyacs.lib.astrotime.decyear2cal(date)
                print(("-- adding XYZ at date %10.5lf mday %02d  month %02d doy %03d " % (date, mday,month,pyacs.lib.astrotime.decyear2dayno(date)[0]) ))
    
    decim_ts.data = keep_data

    # we want the first value to correspond to the first date of the time series
    
    decim_ts.data[0,0] = self.data[0,0]

        
    # forced dates
    
    for forced_date in dates:
        sub_ts = self.extract_periods([1980.0, forced_date])
        
        if sub_ts.data is not None:
            
            decim_ts.data = np.vstack( (decim_ts.data,sub_ts.data[-1,:].reshape(1,-1) ) )
            if verbose:
                new_date = sub_ts.data[-1,0]
                print('-- adding forced date: ', new_date, ' ','-'.join(map(str,pyacs.lib.astrotime.decyear2cal(new_date))[:2]))

    
    # creates .data_xyz
    decim_ts.neu2xyz(corr=True)
    
    # write ts as apr

    # open apr file

    
    if os.path.exists(apr):
        apr_file=open(apr,'a')
    else:
        apr_file=open(apr,'w')

    [cx,cy,cz] = decim_ts.data_xyz[0,1:4]

    for i in np.arange( decim_ts.data.shape[0] ):
        
        current_date = decim_ts.data_xyz[i,0]
        
        (mday,month,_ut)=pyacs.lib.astrotime.decyear2cal( current_date )
        (doy,_ut)=pyacs.lib.astrotime.decyear2dayno( current_date )
        year=int( current_date )
        
        x,y,z = decim_ts.data_xyz[i,1:4]
        
        d = np.sqrt( np.sum( (cx-x)**2 + (cy-y)**2 +(cz-z)**2 ) )
        
        if (dates_gap.shape[0]) > 0 and np.min( np.sqrt( (dates_gap - current_date)**2 ) )  < (1. / 365.25):
            str_comment = ("%s gap  update %d %02d %02d %03d delta %.3lf m" % (str_comment_gen, year,month,mday,doy,d ))
        elif ( dates != [] ) and ( np.min( np.sqrt( (np.array(dates) - current_date)**2) ) < (1. / 365.25)):
            str_comment = ("%s user update %d %02d %02d %03d delta %.3lf m" % (str_comment_gen, year,month,mday,doy,d ))
        else:
            str_comment = ("%s tol  update %d %02d %02d %03d delta %.3lf m" % (str_comment_gen, year,month,mday,doy,d ))
    
        apr_file.write(' %s_GPS %16.4lf  %16.4lf  %16.4lf    0.0000 0.0000 0.0000 %8.4lf %s\n' % \
                       (decim_ts.code.upper(),decim_ts.data_xyz[i,1],decim_ts.data_xyz[i,2],decim_ts.data_xyz[i,3],current_date ,str_comment) )
    
        cx = x
        cy = y
        cz = z
        
    
    apr_file.close()
     


###################################################################
def save_apr(self,apr, epoch=None, verbose=False, excluded_periods=None):
###################################################################
    """
    save results of a Gts analysis in globk format apr file
    
    :param apr: apr file (Globk format)
    :param epoch: epoch in decimal year for coordinates in apr
    :param verbose: verbose mode (boolean)
    :param exluded_periods: periods to be excluded

    :return: Gts instance
    
    :note: following Globk's convention, site will be named XXXX_1PS, XXXX_2PS etc between offset dates
    
    """
    import pyacs.lib.coordinates
    import numpy as np

    import os

    # populates lperiod

    sdate=self.data[0,0]
    edate=self.data[-1,0]
    
    csdate=sdate
    cedate=edate

    lperiod=[]
    
    for offset_date in self.offsets_dates:
        cedate=offset_date
        if verbose:
            print('-- period ',[csdate,cedate],' for site ',self.code)
        lperiod.append([csdate,cedate])
        csdate=cedate
    
    cedate=edate
    if verbose:
        print('-- period ',[csdate,cedate],' for site ',self.code)
    lperiod.append([csdate,cedate])
    
    # open apr file
    
    if os.path.exists(apr):
        apr_file=open(apr,'a')
    else:
        apr_file=open(apr,'w')
    
        
    # loop on periods
    index=0
    
    for period in lperiod:
        
        index = index + 1
        
        # extract time series
        w_gts=self.extract_periods([period])
        
        # calculates the NEU position - use an average of the five first days
        tneu=np.mean(w_gts.data[:4,0:4],axis=0)
        
        # get neu velocity
        try:
            vel=w_gts.detrend().velocity
        except:
            vel=np.zeros(6)
            
        # propagates NEU coordinates at epoch if provided
        if epoch is not None:
            cepoch=epoch
            neu_at_cepoch=tneu[1:]+(epoch-tneu[0])*vel[0:3]
        else:
            cepoch=tneu[0]
            neu_at_cepoch=tneu[1:]
        
        # converts NEU & vel in XYZ
        
        X,Y,Z=pyacs.lib.coordinates.denu_at_x0y0z0_to_xyz(neu_at_cepoch[1],neu_at_cepoch[0],neu_at_cepoch[2],w_gts.X0, w_gts.Y0, w_gts.Z0)
        VX,VY,VZ=pyacs.lib.coordinates.denu_at_x0y0z0_to_xyz(vel[1],vel[0],vel[2],w_gts.X0, w_gts.Y0, w_gts.Z0)

        VX=VX-w_gts.X0
        VY=VY-w_gts.Y0
        VZ=VZ-w_gts.Z0

        
        if index == 1:
            apr_file.write(' %s_GPS %16.4lf  %16.4lf  %16.4lf    0.0000 0.0000 0.0000 %8.4lf %15.5lf  %15.5lf  %15.5lf 0.00000 0.00000 0.00000 from pyacs %s\n' % \
                       (w_gts.code.upper(),X,Y,Z,cepoch,VX,VY,VZ,os.getcwd()))
            
        else:
            apr_file.write(' %s_%1dPS %16.4lf  %16.4lf  %16.4lf    0.0000 0.0000 0.0000 %8.4lf %15.5lf  %15.5lf  %15.5lf 0.00000 0.00000 0.00000 from pyacs %s\n' % \
                       (w_gts.code.upper(),index,X,Y,Z,cepoch,VX,VY,VZ,os.getcwd()))
        
    apr_file.close()


    return(self)

###################################################################
def read_offset_dates(self,offset_file):
###################################################################
    """
    Reads an offset file and populates offsets_dates (pyacs format) attribute of the current Gts instance.
    format is simply a code dates. dates can be any format read by pyacs.guess_date 
    
    :param offset_file: offset_file to be read

    :return: the current Gts instance
    
    """

    import pyacs.lib.astrotime
    
    try:
        fs=open(offset_file,'r')
    except:
        print("-- Could not open ",offset_file) 
        return(None)

    for line in fs:
        if (len(line)<2):continue
        if (line[0]=='#'):continue
        lline=line.split()
        
        if lline[0] == self.code:
            self.offsets_dates.append(pyacs.lib.astrotime.guess_date(lline[-1]))

    return(self)


###################################################################
def info(self,info=2):
###################################################################
    """
    Print various informations about the time series
    """

    import numpy as np
    import pyacs.lib.astrotime
    
    # conversion from p,q to amplitud & phase of seasonal terms

    ###################################################################
    def __pq2aphi__(p,q):
        return(np.sqrt(p**2+q**2), np.arctan2(p,q)/(2.*np.pi))
    

    ###################################################################
    def __print_header__(text):
        print("########################################################################")
        print("# ",text.upper())
        print("########################################################################")


    # period
    time_series_length=self.data[-1,0]-self.data[0,0]
    __print_header__(("general information site: %s" % self.code))
    print(("# Time series period: %9.4lf - %9.4lf (%6.3lf yr)" %(self.data[0,0],self.data[-1,0],time_series_length)))
    sdoy,_ut=pyacs.lib.astrotime.decyear2dayno(self.data[0,0])
    edoy,_ut=pyacs.lib.astrotime.decyear2dayno(self.data[-1,0])
    
    smjd=pyacs.lib.astrotime.decyear2mjd(self.data[0,0])
    emjd=pyacs.lib.astrotime.decyear2mjd(self.data[-1,0])

    print(("# Time series period: %4d  %03d - %4d  %03d  (%d days)" %(int(self.data[0,0]),sdoy,int(self.data[-1,0]),edoy,int(emjd-smjd))))
    
    # obs
    n_obs=self.data[:,0].shape[0]

    n_expected_obs=int(emjd-smjd)+1
    loss=100.0 - float(n_obs)/float(n_expected_obs)*100.0
    print(("# Obs %5d # Expected %5d # Loss %5.1lf %%" %(n_obs,n_expected_obs,loss)))
    
    # approximate coordinate
    if (self.lon,self.lat) != (None,None):
        print(("# App. coordinates (deg.dec): longitude: %8.4lf latitude: %8.4lf " % (self.lon,self.lat)))
    else:
        print ("# App. coordinates (deg.dec): not provided")

    # velocity
    __print_header__("velocity")
    
    if not isinstance(self.velocity,np.ndarray):print ("# Velocity estimates: not available")
    else:
        print ("   V_North     V_East       V_Up     S_V_North    S_V_East      S_V_Up   ")
        if self.velocity.shape[0]==6:
            print(("%10.2lf %10.2lf %10.2lf    %10.2lf  %10.2lf  %10.2lf  " % tuple(self.velocity*1000.0)))
        else:
            print(("%10.2lf %10.2lf %10.2lf        N/A         N/A          N/A " % tuple(self.velocity*1000.0)))

    # annual terms
    __print_header__("seasonal terms")
    if not isinstance(self.annual,np.ndarray):print ("# Annual terms     : not available")
    else:
        print("# Annual and semi-annual terms - phase in years")
        print ("N_amplitude    N_phase  E_amplitude    E_phase   U_amplitude     U_phase")
        a_n,phi_n=__pq2aphi__(self.annual[0],self.annual[1])
        a_e,phi_e=__pq2aphi__(self.annual[2],self.annual[3])
        a_u,phi_u=__pq2aphi__(self.annual[4],self.annual[5])
        
        print((" %10.4lf %10.4lf   %10.4lf %10.4lf    %10.4lf  %10.4lf  (annual)" % (a_n,phi_n,a_e,phi_e,a_u,phi_u)))

    # semi-annual terms
    if not isinstance(self.semi_annual,np.ndarray):print ("# Semi annual terms: not available")
    else:
        a_n,phi_n=__pq2aphi__(self.semi_annual[0],self.semi_annual[1])
        a_e,phi_e=__pq2aphi__(self.semi_annual[2],self.semi_annual[3])
        a_u,phi_u=__pq2aphi__(self.semi_annual[4],self.semi_annual[5])
        
        print((" %10.4lf %10.4lf   %10.4lf %10.4lf    %10.4lf  %10.4lf  (semi-annual)" % (a_n,phi_n,a_e,phi_e,a_u,phi_u)))


    # offsets dates
    __print_header__("offsets")

    if len(self.offsets_dates) == 0:print(("# offsets_dates    :  %4d " % 0))
    else:
        print(("# offsets_dates    :  %4d " % len(self.offsets_dates)))
        if info>0:
            for date in self.offsets_dates:
                doy,_ut= pyacs.lib.astrotime.decyear2dayno(date)
                mday,month,_ut=pyacs.lib.astrotime.decyear2cal(date)
                print(("   %9.4lf # %03d %4d # %02d/%02d/%4d" % (date,doy,int(date),mday,month,int(date))))
    
    # estimated offsets dates
    if not isinstance(self.offsets_values,np.ndarray):
        print(("# estimated offsets:  %4d " % 0))
    else:
        print(("# estimated offsets:        %4d " % self.offsets_values.shape[0]))
        if info>0:
            print(" - Estimated offsets in unit*1000 of original file")

            print("     Date      North        East          Up     S_North      S_East        S_up")
            for index in range(self.offsets_values.shape[0]):
                o_values=self.offsets_values[index,:]*1000.0
                o_values[0]=o_values[0]/1000.0
                print(("%9.4lf %10.2lf  %10.2lf  %10.2lf  %10.2lf  %10.2lf  %10.2lf  " % tuple(o_values)))

    # outliers
    __print_header__("outliers")
    if len(self.outliers) == 0:print(("# outliers:          %4d " % 0))
    else:
        print(("# outliers:         %4d " % len(self.outliers)))
        if info>0: 
            print(" - dates of outliers")
            for index in self.outliers:
                date = self.data[index,0]
                doy,_ut= pyacs.lib.astrotime.decyear2dayno(self.data[index,0])
                mday,month,_ut=pyacs.lib.astrotime.decyear2cal(self.data[index,0])
                print(("   %9.4lf # %03d %4d # %02d/%02d/%4d" % (self.data[index,0],doy,int(date),mday,month,int(date))))
    return(self)
