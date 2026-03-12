"""
Read PRIDE kinematic files (PRIDE-PPPAR format).
"""

###############################################################################
def read_pride(self,tsdir='.',tsfile=None, xyz=True, verbose=False):
###############################################################################
    """
    Read PRIDE-PPPAR kinematic result file.

    Parameters
    ----------
    tsdir : str, optional
        Directory of pride-pppar kinematic files.
    tsfile : str, optional
        Pride-pppar kinematic file to load. If None, looks for kin_*code.
    xyz : bool, optional
        If True, keep data in XYZ; otherwise convert to NEU.
    verbose : bool, optional
        Verbose mode.

    Notes
    -----
    If tsfile is None, read_pride looks for a file named kin_*code.
    """

    # import
    import numpy as np
    import pyacs.lib.astrotime
    from glob import glob
    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG
    import re


    def read_old_format(pride_file):
        """
        Read PRIDE-PPPAR kinematic result file (old format).

        Parameters
        ----------
        pride_file : str
            Path to pride-pppar kinematic file.

        Returns
        -------
        None
            Data are loaded into self.data_xyz and self.data.
        """
        # read pride file
        data = np.genfromtxt(pride_file,skip_header=3)
        # remove lines with star indicating processing problems
        if np.isnan( data ).any():
            WARNING("Nan found in: %s" % pride_file )
            WARNING("Removing lines: " , np.argwhere(np.isnan(data))[:,0].flatten())
            data = data[~np.isnan(data).any(axis=1)]
    
        # fill future .data_xyz
        data_mod = np.zeros((data.shape[0],10))
        data_mod[:,4:7] = 1.E-3
        data_mod[:,0] = data[:,0] + data[:,1] / (60 * 60 * 24. )
    
        data_mod[:,1:4] = data[:,2:]
        if self.data_xyz is None:
            self.data_xyz = data_mod
        else:
            self.data_xyz = np.vstack((self.data_xyz,data_mod))
    
        self.xyz2neu( corr=False , verbose=verbose )
        self.data_xyz[:,0] = self.data[:,0] = pyacs.lib.astrotime.mjd2decyear( self.data_xyz[:,0] )

        # fill t0
        self.t0=self.data[0,0]

        # fill lon, lat
        lon_radian,lat_radian,self.h=pyacs.lib.coordinates.xyz2geo(self.X0,self.Y0,self.Z0)
        self.lon=np.degrees(lon_radian)
        self.lat=np.degrees(lat_radian)
        
        # check duplicate or non-ordered entries
        #if self.data.shape[0]>1:
        #    if np.min(np.diff(self.data[:,0])) <=0:
        #        print("!!! time series not properly ordered by dates or dates duplicated ")
        #        self.reorder()

        # force clean

        self.offsets_dates=[]
        self.offsets_values=None
        self.outliers=[]
        self.annual=None
        self.semi_annual=None
        self.velocity=None

        return(self)


    def read_new_format(filename):
        """Read PRIDE-PPPAR kinematic result file (new format).

        Data are read after the line containing the header
        "* Mjd       Sod               X             Y             Z         Latitude        Longitude        Height".

        Parameters
        ----------
        filename : str
            Path to the PRIDE-PPPAR kinematic file.

        Returns
        -------
        self
            Gts instance with data loaded.
        """
        # pattern
        pattern = '* Mjd       Sod               X             Y             Z         Latitude        Longitude        Height'    
        # First, find the line where the data starts
        with open(filename, 'r') as f:
            lines = f.readlines()

        # Detect the line matching your pattern
        for i, line in enumerate(lines):
            if pattern in line:
                data_start = i + 1  # Assuming data starts on the next line
                break
        else:
            raise ValueError("Pattern not found in file")

        # Read the data part using np.genfromtxt
        data = np.genfromtxt(lines[data_start:], usecols=(0,1,2,3,4), delimiter=None)  # Customize delimiter if needed
        # remove lines with star indicating processing problems
        if np.isnan( data ).any():
            WARNING("Nan found in: %s" % pride_file )
            WARNING("Removing lines: " , np.argwhere(np.isnan(data))[:,0].flatten())
            data = data[~np.isnan(data).any(axis=1)]
    
        # fill future .data_xyz
        data_mod = np.zeros((data.shape[0],10))
        data_mod[:,4:7] = 1.E-3
        data_mod[:,0] = pyacs.lib.astrotime.mjd2decyear( data[:,0] + data[:,1] / (60 * 60 * 24. ) )
        data_mod[:,1:4] = data[:,2:]

        if self.data_xyz is None:
            self.data_xyz = data_mod
        else:
            self.data_xyz = np.vstack((self.data_xyz,data_mod))
    
        self.xyz2neu( corr=False , verbose=verbose )


        # fill t0
        self.t0=self.data[0,0]

        # X0,Y0,Z0
        self.X0,self.Y0,self.Z0 = self.data_xyz[0,2:5]

        # fill lon, lat
        lon_radian,lat_radian,self.h=pyacs.lib.coordinates.xyz2geo(self.X0,self.Y0,self.Z0)
        self.lon=np.degrees(lon_radian)
        self.lat=np.degrees(lat_radian)
        
        # check duplicate or non-ordered entries
        #if self.data.shape[0]>1:
        #    if np.min(np.diff(self.data[:,0])) <=0:
        #        print("!!! time series not properly ordered by dates or dates duplicated ")
        #        self.reorder()

        # force clean

        self.offsets_dates=[]
        self.offsets_values=None
        self.outliers=[]
        self.annual=None
        self.semi_annual=None
        self.velocity=None

        return(self)


    # name of the file to be read - if not provided, tries to guess

    # READ DIR OF USER SPECIFIED FILE     
    if (tsfile is None):
        if (self.code is not None):
            try:
                l_pride_file=glob( tsdir+'/kin_*'+self.code.lower() )
            except:
                ERROR(f"Could not find any time series file for code {self.code}")
                return()
        else:
            ERROR("no code or file provided.")

    else:
        l_pride_file = [ tsfile]

    for pride_file in sorted( l_pride_file ):
        VERBOSE("reading: %s " % pride_file)
        try:
            read_new_format(pride_file)
        except:
            try:
                read_old_format(pride_file)
            except:
                ERROR(f"Error reading {pride_file}. Returning None.")
                return(None)



###############################################################################
def read_pride_pos(self,tsdir='.',tsfile=None, verbose=False):
###############################################################################
    """Read PRIDE-PPPAR static result file(s).

    Parameters
    ----------
    tsdir : str, optional
        Directory containing pride-pppar pos static files. Default '.'.
    tsfile : str, optional
        Specific file to load. If None, searches for pos_*code in tsdir.
    verbose : bool, optional
        If True, enable verbose output.

    Notes
    -----
    If tsfile is None, looks for files named pos_*<code> in tsdir.
    """

    # import
    import numpy as np
    import pyacs.lib.astrotime
    from glob import glob
    import pyacs.lib.glinalg
    
    # name of the file to be read - if not provided, tries to guess
    
    if (tsfile is None):
        if (self.code is not None):

            try:
                l_pride_file=glob( tsdir+'/pos_*'+self.code.lower() )
            except:
                print("!!! Error: Could not find any time series file for code ",self.code)
                return()
        else:
            print("!!! Error: no code or file provided.")

    else:
        l_pride_file = [ tsfile]

    if verbose:
        print('-- will read: ')
        for pride_file in sorted( l_pride_file ):
            print("%s" % pride_file)
    for pride_file in sorted( l_pride_file ):
        if verbose:
            print("-- reading: %s " % pride_file)
        
        # read pride file
        d=np.genfromtxt(pride_file,skip_header=2,skip_footer=1)              
        [X,Y,Z] = d[0,1:4]
        [SX,SY,SZ] = sigma_m = d[1,1:4]
        corr_coef = d[1,1:4]
        CORR = np.eye((3))
        CORR[0,1] = CORR[1,0] = corr_coef[0]
        CORR[0,2] = CORR[2,0] = corr_coef[1]
        CORR[1,2] = CORR[2,1] = corr_coef[2]

        COV = pyacs.lib.glinalg.corr_to_cov(CORR, sigma_m)
        # get the date
        fp = open(pride_file)
        mjd = float( fp.readline().split()[-1])
        decyear =  pyacs.lib.astrotime.mjd2decyear( mjd )
        # array to stak
        
        obs_array = np.array([[decyear,X,Y,Z,SY,SY,SZ,COV[0,1],COV[0,2],COV[1,2]]])
        
        # fill  .data_xyz
      
        if self.data_xyz is None:
            self.data_xyz = obs_array
        else:
            self.data_xyz = np.vstack((self.data_xyz,obs_array))
    
    self.xyz2neu( corr=False , verbose=verbose )

    # fill t0
    self.t0=self.data[0,0]

    # fill lon, lat
    lon_radian,lat_radian,self.h=pyacs.lib.coordinates.xyz2geo(self.X0,self.Y0,self.Z0)
    self.lon=np.degrees(lon_radian)
    self.lat=np.degrees(lat_radian)

    self.offsets_dates=[]
    self.offsets_values=None
    self.outliers=[]
    self.annual=None
    self.semi_annual=None
    self.velocity=None

    return(self)
