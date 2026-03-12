"""Read and write GAMIT/GLOBK PBO pos files."""


###################################################################
## Loads Gts from GAMIT/GLOBK  pos file format for time series
###################################################################
def read_pos(self,tsdir='.',tsfile=None, xyz=True, verbose=False):
    """Read GAMIT/GLOBK PBO pos file and load the time series into this Gts.

    Parameters
    ----------
    tsdir : str, optional
        Directory containing pos file(s). Default is '.'.
    tsfile : str, optional
        Pos file to load. If None, a file CODE*.pos is sought. Default is None.
    xyz : bool, optional
        If True, read XYZ and sx,sy,sz, corr columns. Default is True.
    verbose : bool, optional
        If True, print progress. Default is False.

    Returns
    -------
    Gts
        self (data, code, X0, Y0, Z0, t0 populated).

    Notes
    -----
    A pos file contains (almost) all needed info. If tsfile is None,
    read_pos looks for a file named CODE*.pos.
    """

    # import
    import numpy as np
    import pyacs.lib.astrotime
    import pyacs.lib.coordinates
    import os


    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG

    from icecream import ic

    # name of the file to be read - if not provided, tries to guess
    
    if (tsfile is None):
        if (self.code is not None):
            from glob import glob
            try:
                pos_file=glob(tsdir+'/*'+self.code.upper()+'*.pos')[0]
            except:
                ERROR("Could not find any time series file for code ",self.code, exit=True)
        else:
            ERROR("no code or file provided.", exit=True)

    else:
        pos_file=tsfile

    # file

    self.ifile=os.path.abspath(pos_file)

    # actual read

    # get site code from header, handling both 4-character and 9-character IDs
    import re
    self.code = None
    with open(pos_file, "r") as f_pos:
        for _ in range(40):
            line = f_pos.readline()
            if not line:
                break
            match = re.match(r"\s*(?:4|9)-character ID\s*:\s*(\S+)", line, flags=re.IGNORECASE)
            if match:
                header_id = match.group(1).strip()
                self.code = header_id[:4].upper() if len(header_id) >= 4 else header_id.upper()
                break
    if self.code is None:
        # fallback for non-standard headers: keep backward compatibility with former behavior
        import linecache
        header_id_line = linecache.getline(pos_file, 3)
        header_id = header_id_line.split(':', 1)[-1].strip() if ':' in header_id_line else header_id_line.strip()
        self.code = header_id[:4].upper() if len(header_id) >= 4 else header_id.upper()

    # determine the header length
    # last header line has *YYYYMMDD

    header_length = 0
    with open(pos_file, "r") as f_pos:
        for line in f_pos:
            if line.strip() and not line.startswith('*YYYYMMDD'):
                header_length += 1
            else:
                header_length += 1
                break
    if header_length == 0:
        ERROR("Could not determine header length. Returning Gts", exit=True)


    if xyz:
    # trust xyz and re-creates dneu
        data=np.genfromtxt(pos_file,skip_header=header_length,usecols=tuple(range(12)))
        # reshape to ensure a 2D array
        if data.ndim == 1:
            data=data.reshape((1,data.shape[0]))

        # convert dates to dec.year
        mjd_data=list(data[:,2])
        dec_year=list(map(pyacs.lib.astrotime.mjd2decyear,mjd_data))

        # fill data_xyz
        self.data_xyz=data[:,2:12]
        
        # fill data
        self.data_xyz[:,0]=np.array(dec_year)
        self.xyz2neu(corr=True)
        
        if np.isnan(self.data[:,4:]).any():
            self.data[:,4:7]=1.E-3
            self.data[:,7:10]=0.
        
    else:
    # for some reason, you prefer to trust dneu
        try:
            data=np.genfromtxt(pos_file,skip_header=header_length,usecols=tuple(range(15,24)))
            # reshape to ensure a 2D array
            if data.ndim==1:
                data=data.reshape((1,data.shape[0]))
        except:
            # there is a big problem
            raise IOError("!!! Error: Could not read dN,dE,dU from file: %s " % pos_file )
        
        # convert dates to dec.year
        mjd_data=list(np.genfromtxt(pos_file,skip_header=header_length,usecols= 2))
        dec_year=pyacs.lib.astrotime.mjd2decyear(mjd_data)

        # XYZ reference from header
        import linecache
        [self.X0,self.Y0,self.Z0]=list(map(float , linecache.getline(pos_file, 8).split(':')[-1].split()[:3] ))

        # fill data
        
        self.data = np.zeros( (dec_year.shape[0],10) )
        
        self.data[:,1:]=data
        self.data[:,0]=np.array(dec_year)
    
    # fill t0
    self.t0=self.data[0,0]

    # fill lon, lat
    lon_radian,lat_radian,self.h=pyacs.lib.coordinates.xyz2geo(self.X0,self.Y0,self.Z0)
    self.lon=np.degrees(lon_radian)
    self.lat=np.degrees(lat_radian)
    
    # check duplicate or non-ordered entries
    if self.data.shape[0]>1:
        if np.min(np.diff(self.data[:,0])) <=0:
            print("!!! time series not properly ordered by dates or dates duplicated ")
            self.reorder()

    # force clean

    self.offsets_dates=[]
    self.offsets_values=None
    self.outliers=[]
    self.annual=None
    self.semi_annual=None
    self.velocity=None

    return(self)


###################################################################
def write_pos(self,idir='./pos',add_key='' , force=None, verbose=False):
###################################################################
    """Write time series in GAMIT/GLOBK PBO pos format.

    Parameters
    ----------
    idir : str, optional
        Output directory. Default is './pos'.
    add_key : str, optional
        If non-blank, output file is CODE_add_key.pos; otherwise CODE.pos.
    force : str, optional
        'data' or 'data_xyz' to force source; None for default behavior.
    verbose : bool, optional
        If True, print progress. Default is False.

    Returns
    -------
    Gts
        self.

    Notes
    -----
    Default (force=None): if both data and data_xyz exist, both are written;
    if only data, uses X0,Y0,Z0 to write data_xyz; if only data_xyz, recreates
    data and writes.
    """

    # import
    import pyacs.lib.coordinates
    import pyacs.lib.astrotime
    import numpy as np
    import os

    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG


    # force option
    
    if force is not None:
        if force not in ['data','data_xyz']:
            ERROR('force parameter must be either data or data_xyz. Returning Gts')
            return(self)
    
        if force == 'data':
            self.data_xyz = None

        if force == 'data_xyz':
            self.data = None

    # cdata
    
    if force != 'data_xyz':
        if not self.cdata(data=True):
            WARNING('Can not write .pos file. Problem with .data attribute. Returning Gts')
            self.cdata(data=True,verbose=True)
            return(self)


    ###############################################
    def __print_header_pos_file(f_pos,code,first_epoch,last_epoch,release_date,XYZ_ref,geo_ref,verbose=False):
        f_pos.write("PBO Station Position Time Series. Reference Frame : Unknown\n")
        f_pos.write("Format Version: 1.1.0\n")
        f_pos.write("4-character ID: %s\n" % code)
        f_pos.write("Station name  : %s        \n" % code)
        f_pos.write("First Epoch   : %s\n" % first_epoch)
        f_pos.write("Last Epoch    : %s\n" % last_epoch)
        f_pos.write("Release Date  : %s\n" % release_date)
        f_pos.write("XYZ Reference position :   %s\n" % XYZ_ref)
        f_pos.write("NEU Reference position :   %s (Unknown/WGS84)\n" % geo_ref)
        f_pos.write("Start Field Description\n")
        f_pos.write("YYYYMMDD      Year, month, day for the given position epoch\n")
        f_pos.write("HHMMSS        Hour, minute, second for the given position epoch\n")
        f_pos.write("JJJJJ.JJJJJ   Modified Julian day for the given position epoch\n")
        f_pos.write("X             X coordinate, Specified Reference Frame, meters\n")
        f_pos.write("Y             Y coordinate, Specified Reference Frame, meters\n")
        f_pos.write("Z             Z coordinate, Specified Reference Frame, meters\n")
        f_pos.write("Sx            Standard deviation of the X position, meters\n")
        f_pos.write("Sy            Standard deviation of the Y position, meters\n")
        f_pos.write("Sz            Standard deviation of the Z position, meters\n")
        f_pos.write("Rxy           Correlation of the X and Y position\n")
        f_pos.write("Rxz           Correlation of the X and Z position\n")
        f_pos.write("Ryz           Correlation of the Y and Z position\n")
        f_pos.write("Nlat          North latitude, WGS-84 ellipsoid, decimal degrees\n")
        f_pos.write("Elong         East longitude, WGS-84 ellipsoid, decimal degrees\n")
        f_pos.write("Height (Up)   Height relative to WGS-84 ellipsoid, m\n")
        f_pos.write("dN            Difference in North component from NEU reference position, meters\n")
        f_pos.write("dE            Difference in East component from NEU reference position, meters\n")
        f_pos.write("du            Difference in vertical component from NEU reference position, meters\n")
        f_pos.write("Sn            Standard deviation of dN, meters\n")
        f_pos.write("Se            Standard deviation of dE, meters\n")
        f_pos.write("Su            Standard deviation of dU, meters\n")
        f_pos.write("Rne           Correlation of dN and dE\n")
        f_pos.write("Rnu           Correlation of dN and dU\n")
        f_pos.write("Reu           Correlation of dEand dU\n")
        f_pos.write("Soln          corresponding to products  generated with rapid or final orbit products, in supplemental processing, campaign data processing or reprocessing\n")
        f_pos.write("End Field Description\n")
        f_pos.write("*YYYYMMDD HHMMSS JJJJJ.JJJJ         X             Y             Z            Sx        Sy       Sz     Rxy   Rxz    Ryz            NLat         Elong         Height         dN        dE        dU         Sn       Se       Su      Rne    Rnu    Reu  Soln\n")

    ###############################################
    def __dndedu_to_pos_file(gts,verbose=False):
        """
        Assuming that the Gts instance has dn de du coordinates, returns an list ready to be written a Gamit/Globk .pos file
        Requires the Gts instance to have high accuracy .lon .lat & .h
        Since X Y Z sigmas cannot be properly calculated they are set to 0.0
        """

        #################################
        #def __format_cal_pos(lyear,lcal):
        #    lformatted_date1=[]
        #    lformatted_date2=[]
        #    
        #    for i in range(len(lyear)):
        #        lformatted_date1.append(int("%d%02d%02d" % (lyear[i],lcal[i][1],lcal[i][0])))
        #        uts=pyacs.lib.astrotime.ut2uts(lcal[i][2])
        #        (h,mn,s,microsecond)=pyacs.lib.astrotime.uts2hmsmicros(np.around(uts))
        #        lformatted_date2.append(int(("%02d%02d%02d" % (h,mn,s))))
        #    return(lformatted_date1,lformatted_date2)

        #################################
        def __format_cal_pos(np_decyear):
            lformatted_date1=[]
            lformatted_date2=[]
            
            for i in range(len(np_decyear)):
                
                mdt = pyacs.lib.astrotime.datetime_round_second( pyacs.lib.astrotime.decyear2datetime(np_decyear[i]) )
                
                lformatted_date1.append(int("%d%02d%02d" % (mdt.year,mdt.month,mdt.day)))
                lformatted_date2.append(int(("%02d%02d%02d" % (mdt.hour,mdt.minute,mdt.second))))
            return(lformatted_date1,lformatted_date2)

        
        # import
        import pyacs.lib.coordinates
        import pyacs.lib.astrotime
        import numpy as np
        

        # check for duplicate records
        VERBOSE('running reorder')
        self.reorder()

        dndedulamphih=np.zeros((gts.data.shape[0],6))
        dndedulamphih[:,:3]=self.data[:,1:4]

        # case dndedu and XYZ exists
        if isinstance(gts.data,np.ndarray) and isinstance(gts.data_xyz,np.ndarray):
            VERBOSE("will write NEU data & XYZ data independently")
            # write data & data_xyz independently
            XYZ=gts.data_xyz[:,1:10]
            
        # case dndedu exits and XYZ does not exists
        if isinstance(gts.data,np.ndarray) and not isinstance(gts.data_xyz,np.ndarray):
            # creates data_xyz independently
            VERBOSE('creating XYZ data from NEU')
            
            if gts.X0 is not None:
                dndedulamphih[:,3]=gts.X0
                dndedulamphih[:,4]=gts.Y0
                dndedulamphih[:,5]=gts.Z0
                # CHECK ENU or NEU
                tmp_XYZ=list(map(pyacs.lib.coordinates.denu_at_x0y0z0_to_xyz,dndedulamphih[:,1],dndedulamphih[:,0],dndedulamphih[:,2],dndedulamphih[:,3],dndedulamphih[:,4],dndedulamphih[:,5]))
                tmp_XYZ=np.array(tmp_XYZ)
                XYZ=np.zeros((gts.data.shape[0],10))
                XYZ[:,:3]=tmp_XYZ[:,:3]
                
            else:
                # put 0 values
                XYZ=np.zeros((gts.data.shape[0],10))

        # case dndedu does exist and XYZ exists
        if isinstance(gts.data_xyz,np.ndarray) and not isinstance(gts.data,np.ndarray):
            # write data & data_xyz independently
            VERBOSE('creating NEU data from XYZ')
            gts.xyz2neu(corr=True)
            XYZ=gts.data_xyz[:,1:10]
        
        # first & second columns - calendar date and time
        #lcal=map(list,map(pyacs.lib.astrotime.decyear2cal,gts.data[:,0]))
        #lyear=map(int,gts.data[:,0])
        #ldate1,ldate2=__format_cal_pos(lyear,lcal)
        
        ldate1,ldate2=__format_cal_pos(gts.data[:,0])
        
        
        # 3rd column - modified julian day
        lmjd=list(map(pyacs.lib.astrotime.decyear2mjd,gts.data[:,0]))
                
        # column 13-15 - geographical coordinates
        l_13_15=list(map(list,list(map(pyacs.lib.coordinates.xyz2geo,XYZ[:,0],XYZ[:,1],XYZ[:,2]))))

        # other columns        
        lpos=[]
        for i in range(gts.data.shape[0]):
            if self.data.shape[1]==10:
            # self.data has date, dn,de,du,sdn,sde,sdu,Rne,Rnu,Reu
                new_pos=[ldate1[i],ldate2[i],lmjd[i]]+[XYZ[i,0],XYZ[i,1],XYZ[i,2],XYZ[i,3],XYZ[i,4],XYZ[i,5],XYZ[i,6],XYZ[i,7],XYZ[i,8]]+l_13_15[i]+gts.data[i,1:10].tolist()+['pyacs']
            else:
            # self.data has date, dn,de,du,sdn,sde,sdu
                new_pos=[ldate1[i],ldate2[i],lmjd[i]]+[XYZ[i,0],XYZ[i,1],XYZ[i,2],XYZ[i,3],XYZ[i,4],XYZ[i,5],XYZ[i,6],XYZ[i,7],XYZ[i,8]]+l_13_15[i]+gts.data[i,1:7].tolist()+[0.0,0.0,0.0,'pyacs']

            lpos.append(new_pos)
        
        return(lpos)

    
    ###############################################
    def __print_content_pos_file(f_pos,lformatted):

        import numpy as np

        for i in range(len(lformatted)):
            
            [YYYYMMDD,HHMMSS,MJD,X,Y,Z,Sx,Sy,Sz,Rxy,Rxz,Ryz,Elong,NLat,Height,dN,dE,dU,Sn,Se,Su,Rne,Rnu,Reu,Soln]=lformatted[i]
            f_pos.write\
            (" %8d %6d %10.4lf %14.5lf %14.5lf %14.5lf %8.5lf %8.5lf %8.5lf %6.3lf %6.3lf %6.3lf     %14.10lf  %14.10lf %10.5lf  %10.5lf%10.5lf%10.5lf   %8.5lf %8.5lf %8.5lf %6.3lf %6.3lf %6.3lf %s\n" % \
            (YYYYMMDD,HHMMSS,MJD,X,Y,Z,Sx,Sy,Sz,Rxy,Rxz,Ryz,np.degrees(NLat),np.degrees(Elong),Height,dN,dE,dU,Sn,Se,Su,Rne,Rnu,Reu,Soln))

    ###############################################


    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG


    if not os.path.isdir(idir):
        try:
            os.makedirs(idir, exist_ok=True)
        except:
            WARNING("Could not create directory ",idir)
    
    # add_key
    if add_key != '':
        f_pos_name=idir+'/'+self.code+'_'+add_key+'.pos'
    else:
        f_pos_name=idir+'/'+self.code+'.pos'
        
    # open pos file

    if verbose:
        print("-- writing %s" % os.path.abspath(f_pos_name))

    f_pos=open(f_pos_name,'w+')
    
    # X0,Y0,Z0
    if self.X0 is not None:
        pass
    else:
        if self.lon is not None:
            self.X0,self.Y0,self.Z0=pyacs.lib.coordinates.geo2xyz(np.radians(self.lon),np.radians(self.lat),self.h)
        else:
            print("!!! Cannot produce a proper pos file because a reference position is needed")
            print("!!! Setting all XYZ values to 0")

    # populates info for pos file
    
    wpos=__dndedu_to_pos_file(self,verbose=verbose)
            
            
    sdecyear=self.data[0,0]
    (mday,month,ut)=pyacs.lib.astrotime.decyear2cal(sdecyear)
    uts=pyacs.lib.astrotime.ut2uts(ut)
    (h,mn,s, _microsecond)=pyacs.lib.astrotime.uts2hmsmicros(uts)
    first_epoch=("%d%02d%02d %02d%02d%02d"% (int(sdecyear),month,mday,h,mn,s))

    edecyear=self.data[0,0]
    (mday,month,ut)=pyacs.lib.astrotime.decyear2cal(edecyear)
    uts=pyacs.lib.astrotime.ut2uts(ut)
    (h,mn,s, _microsecond)=pyacs.lib.astrotime.uts2hmsmicros(uts)
    last_epoch=("%d%02d%02d %02d%02d%02d"% (int(edecyear),month,mday,h,mn,s))

    import datetime
    
    current_date=datetime.datetime.today()
    
    release_date=("%d%02d%02d %02d%02d%02d"% (current_date.year,current_date.month,current_date.day,current_date.hour,current_date.minute,current_date.second))
    
    if self.X0 is not None:
        XYZ_ref=(" %14.5lf %14.5lf %14.5lf (unkwnon)" % (self.X0,self.Y0,self.Z0))
    else:
        XYZ_ref='unknown'

    if self.lon is not None:        
        geo_ref=(" %16.10lf %16.10lf %16.10lf (unkwnon)" % ( self.lat,self.lon,self.h ))
    else:
        geo_ref='unknown'
    
    __print_header_pos_file(f_pos,self.code,first_epoch,last_epoch,release_date,XYZ_ref,geo_ref)
    __print_content_pos_file(f_pos,wpos)
    
    f_pos.close()
    return(self)

