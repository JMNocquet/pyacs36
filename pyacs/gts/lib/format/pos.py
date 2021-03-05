"""
Reads and write PBO pos files
"""


###################################################################
## Loads Gts from GAMIT/GLOBK  pos file format for time series
###################################################################
def read_pos(self,tsdir='.',tsfile=None, xyz=True, verbose=False):
    """
    Read GAMIT/GLOBK PBO pos file in a directory and actually loads the time series

    :param tsdir: directory of pos file
    :param tsfile: pos file to be loaded
    :param xyz: reads xyz sx sy sz corr_xy corr_xz corr_yz columns
    :param verbose: verbose mode
    :note: Since a pos file includes (almost) all the information, data, code, X0,Y0,Z0,t0 will be populated
    :note: If tsfile=None, then read_pos will look for a file named CODE*.pos
    """

    # import
    import numpy as np
    import pyacs.lib.astrotime
    import pyacs.lib.coordinates
    import os



    # name of the file to be read - if not provided, tries to guess
    
    if (tsfile is None):
        if (self.code is not None):
            from glob import glob
            try:
                pos_file=glob(tsdir+'/'+self.code.upper()+'*.pos')[0]
            except:
                print("!!! Error: Could not find any time series file for code ",self.code)
                return()
        else:
            print("!!! Error: no code or file provided.")

    else:
        pos_file=tsfile

    # file

    self.ifile=os.path.abspath(pos_file)

    # actual read

    # get site code
    import linecache
    self.code = linecache.getline(pos_file, 3).split(':')[-1].strip()

    
    if xyz:
    # trust xyz and re-creates dneu
        data=np.genfromtxt(pos_file,skip_header=37,usecols=tuple(range(12)))
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
            data=np.genfromtxt(pos_file,skip_header=37,usecols=tuple(range(15,24)))
            # reshape to ensure a 2D array
            if data.ndim==1:
                data=data.reshape((1,data.shape[0]))
        except:
            # there is a big problem
            raise IOError("!!! Error: Could not read dN,dE,dU from file: %s " % pos_file )
        
        # convert dates to dec.year
        mjd_data=list(np.genfromtxt(pos_file,skip_header=37,usecols= 2))
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
def write_pos(self,idir,add_key='' , force=None, verbose=False):
###################################################################
    """
    Write a time series in GAMIT/GLOBK PBO pos format

    :param idir: output directory
    :param add_key: if not blank then the output pos file will be CODE_add_key.pos, CODE.pos otherwise.
    :param force: set force to 'data' or 'data_xyz' to force pos to be written from .data or .data_xyz

    :note1:default behaviour (force = None)
        if data and data_xyz are not None, then print them independently
        if there are data only, then uses X0,Y0,Z0 to write data_xyz
        if there are data_xyz only, recreate data and write it
    
    """

    # force option
    
    if force is not None:
        if force not in ['data','data_xyz']:
            print('!!! ERROR: force keyword must be either data or data_xyz')
            return(self)
    
        if force == 'data':
            self.data_xyz = None

        if force == 'data_xyz':
            self.data = None

    # cdata
    
    if force != 'data_xyz':
        if not self.cdata(data=True):
            print('! Can not write .pos file. Problem with .data attribute.')
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
        
#        # data
#        data=np.copy(gts.data)
#        Hdata = dict(zip(gts.data[:,0],gts.data[:,1:]))
        
        # check for duplicate records
        self.reorder()
        
#        if len(Hdata.keys()) != data.shape[0]:
#            print "!!! dates probably duplicated for site ",gts.code
#            data=np.zeros((len(Hdata.keys()),data.shape[1]))
#        
#        # 
#        data[:,0]=sorted(Hdata.keys())
#        data[:,1:]=collections.OrderedDict(sorted(Hdata.items())).values()
#        
#        gts.data=data
        
        
#        pos=np.zeros((gts.data.shape[0],25))

        dndedulamphih=np.zeros((gts.data.shape[0],6))
        dndedulamphih[:,:3]=self.data[:,1:4]

        # case dndedu and XYZ exists
        if isinstance(gts.data,np.ndarray) and isinstance(gts.data_xyz,np.ndarray):
            if verbose:print("-- write NEU data & XYZ data independently")
            # write data & data_xyz independently
            XYZ=gts.data_xyz[:,1:10]
            
        # case dndedu exits and XYZ does not exists
        if isinstance(gts.data,np.ndarray) and not isinstance(gts.data_xyz,np.ndarray):
            # creates data_xyz independently
            if verbose:print('-- create XYZ data from NEU')
            
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
            if verbose:print('-- create NEU data from XYZ')
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

    import pyacs.lib.coordinates
    import pyacs.lib.astrotime
    import numpy as np
    import os

    # output directory
    if not os.path.isdir(idir):
        try:
            os.mkdir(idir)
            print("-- Creating directory ",idir)
        except:
            print("!!! Error : Could not create directory ",idir)
    
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

