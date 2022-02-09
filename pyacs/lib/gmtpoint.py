
class GMT_Point:
    
    
    """
    Simple Point Class
    
    Attributes (mandatory)
    lon: longitude (decimal degrees)
    lat: latitude (decimal degrees)
    
    Attributes (optinal):
    code : 4-letters code
    he : height above the ellipsoid
    ve, vn: east & north components of velocity 
    sve,svn,sven: formal error (standard deviation) on velocity component and correlation [-1,1] """


###################################################################
    def __init__ (self, code=None,lon=None,lat=None,he=0, Ve=None,Vn=None,Vu=0,\
                  SVe=0,SVn=0,SVu=0,SVen=0,\
                  Cv_xyz=None, Cv_enu=None, index=None):
###################################################################
        try:
            self.code=code
            self.Ve=Ve
            self.Vn=Vn
            self.Vu=Vu
            self.SVe=SVe
            self.SVn=SVn
            self.SVen=SVen
            self.SVu=SVu
            self.lon=lon
            self.lat=lat
            self.he=he
            self.Cv_xyz=Cv_xyz
            self.Cv_enu=Cv_enu
            self.index=index
        except:
            raise ValueError('!!! Error. Could not init GMT_Point.')
        
        
###################################################################
    def copy(self):
###################################################################
        NEW=GMT_Point()
        NEW.code=self.code
        NEW.Ve=self.Ve
        NEW.Vn=self.Vn
        NEW.Vu=self.Vu
        NEW.SVe=self.SVe
        NEW.SVn=self.SVn
        NEW.SVen=self.SVen
        NEW.SVu=self.SVu
        NEW.lon=self.lon
        NEW.lat=self.lat
        NEW.he=self.he
        NEW.Cv_xyz=self.Cv_xyz
        NEW.Cv_enu=self.Cv_enu
        NEW.index=self.index
        
        return(NEW)
        
        
###################################################################
    def get_info(self, display=False, legend=False, verbose=False):
###################################################################
        """
        Print information of the current GMT_Point
        """

        info_legend="code     "
        info=("%4s " % self.code)
        
        info_legend=info_legend + " long.        lat."
        info_legend=info_legend + "       Ve.        Vn."
        info_legend=info_legend+"      SVe.       SVn.       SVen"
        
        info_legend=info_legend+"      V_magnitud    V_azimuth"

        info =info + ("%10.5lf  %10.5lf" % (self.lon,self.lat))

        if ( self.Ve is not None):
            info=info + ("%10.3lf %10.3lf" % (self.Ve,self.Vn))
            (mag,az)=self.magaz()

            if ( self.SVe is not None):
                info=info + ("%10.3lf %10.3lf %10.3lf" % (self.SVe,self.SVn,self.SVen))
                info=info + ("      %10.3lf   %10.1lf " % (mag,az))

        if legend:
            print(info_legend)
            print(info)
        
        if display:
            print(info)
            
        if verbose:
            print('-- returning info string for site: %s' , self.code )
        
        return(info)

###################################################################
    def magaz(self):
###################################################################
        """
        Returns norm & azimuth of velocity from the current GMT_Point
        
        :return: norm,azimuth in mm/yr and decimal degree clockwise from North
        
        :note: sigmas not implemented yet
        
        """
        import math
        # magnitude & azimuth
        mag=math.sqrt(self.Ve**2+self.Vn**2)
#!!!!!!!!!!!!!!!!! NOT CHECKED
        #if (self.SVe != None):
        az=math.degrees(math.atan2(self.Ve,self.Vn))
        return(mag,az)
    
        
###################################################################
    def assign_index(self,index):
###################################################################
        """
        Assigns an index to the current GMT_Point.
        Useful for building linear systems.
        
        :param index: index
        
        :return : the modified GMT_Point instance
        """
        
        self.index=index
        return(self)


###################################################################
    def get_index(self):
###################################################################
        """
        Returns index of the current GMT_Point
        
        :return: index
        """
        return self.index

###################################################################
    def spherical_distance(self,M):
###################################################################

        """
        Returns spherical distance between two GMT_Point (in meters)
        
        :param M: GMT_Point instance
        
        :return : distance along the sphere between current instance and M
        
        """

        import pyacs.lib.coordinates
        import pyacs.message.message as MESSAGE
        import pyacs.message.verbose_message as VERBOSE
        import pyacs.message.error as ERROR
        import pyacs.message.warning as WARNING
        import pyacs.message.debug_message as DEBUG

        DEBUG( ("%.5lf %.5lf %.5lf %.5lf ") % (self.lon, self.lat,  M.lon, M.lat) )

        return pyacs.lib.coordinates.geo_spherical_distance(
            self.lon, self.lat, 0.0, M.lon, M.lat, 0.0, unit='dec_deg')

    
            
###################################################################
    def pole(self,W=None,SW=None,type_euler='rot',option='predict'):
###################################################################

        """Substract the prediction of a given Euler pole"""
        
        import pyacs.lib.coordinates as Coordinates
        import pyacs.lib.euler
        import numpy as np
        
        
        ltype=['euler','rot'] # type is either Euler pole (long. lat in dec. degrees and deg/Myr or cartesian rotation rate vector in rad/yr)
        if (type_euler not in ltype):
            raise ValueError(
                "!!! Error. type_euler should be either euler " +
                "or rot (Euler pole (long. lat in dec. degrees and deg/Myr)" +
                " or cartesian rotation rate vector in rad/yr. type_euler = %s )"
                , type_euler )

        loption=['remove','predict','add'] #
        if (option not in loption):
            raise ValueError(
            "!!! Error. option should be either predict, remove or  add. option = %s " ,
            option
            )
            
        
        if (type_euler == 'euler'):
            if (len(W.shape)!=2):
                W=W.reshape(3,1)

            # convert to rotation rate vector
            (lon, lat, omega)=(W[0,0],W[1,0],W[2,0])
            (wx,wy,wz)=pyacs.lib.euler.euler2rot(lon, lat, omega)
            
            # W
            W=np.zeros((3,1))
            W[0,0]=wx
            W[1,0]=wy
            W[2,0]=wz


        # get spherical coordinates
        
        (x,y,z)=Coordinates.geo2xyz(self.lon,self.lat,self.he,unit='dec_deg')
        (l,p,_)=Coordinates.xyz2geospheric(x,y,z,unit='radians')
        R=Coordinates.mat_rot_general_to_local(l,p)

        # Observation equation in local frame
                
        Ai=np.zeros([3,3],float)
        Ai[0,1]= z
        Ai[0,2]=-y
        Ai[1,0]=-z
        Ai[1,2]= x
        Ai[2,0]= y
        Ai[2,1]=-x
                
        RAi=np.dot(R,Ai)
            
        # Predicted velocity
                
        Pi=np.dot(RAi,W)
            
        pve=Pi[0,0]*1.E3
        pvn=Pi[1,0]*1.E3
        
        N=self.copy()
            
        if (option=='remove'):                
            N.Ve=self.Ve-pve
            N.Vn=self.Vn-pvn
        if (option=='add'):                
            N.Ve=self.Ve+pve
            N.Vn=self.Vn+pvn
        if (option=='predict'):                
            N.Ve=pve
            N.Vn=pvn

        if SW != None:
            # uses pole variance-covariance matrix
            VCV=SW
            #print VCV.shape,Ai.shape
#            VCV=np.zeros((3,3))
#            VCV[0,0]=SW[0]
#            VCV[1,1]=SW[1]
#            VCV[2,2]=SW[2]
        
#            VCV[1,2]=SW[3]
#            VCV[1,3]=SW[4]
#            VCV[2,3]=SW[5]
            
#            VCV_PREDICTION=np.dot(np.dot(Ai,VCV),Ai.T)
            VCV_PREDICTION_ENU=np.dot(np.dot(RAi,VCV),RAi.T)
            (N.SVe,N.SVn)=(np.sqrt(VCV_PREDICTION_ENU[0,0])*1000.0,np.sqrt(VCV_PREDICTION_ENU[1,1])*1000.0 )  

        return(N)
        
###################################################################
    def add_to_gmt_psvelo(self,fname,overwrite=False,verbose=False):
###################################################################

        """
        Adds the current GMT_Point to a gmt psvelo file
        
        :param fname: gmt psvelo file
        :param overwrite: will overwrite the line corresponding having the same code (default is False), generate an error otherwises
        :param verbose: verbose mode
        """

        # reads the velocity file
        import pyacs.lib.vel_field
        
        vel = pyacs.lib.vel_field.Velocity_Field()
        try:
            vel.read(file_name=fname,verbose=verbose)
        except:
            if verbose:
                print('-- file does not exist. Creating: %s ' , fname )

        # if overwrite option
        
        if overwrite:
            if self.code in vel.lcode():
                if verbose:
                    print('-- replacing site %s:%s' % (self.code,vel.site(self.code).get_info()))
                    print('-- with      site %s:%s' % (self.code,self.get_info()))

            vel.remove(self.code)

        # add self
        
        vel.add_point(self)

        # write output psvelo file
        
        vel.write(fname, verbose=verbose)

###################################################################
    def rotate_vel(self, angle,unit='radians'):
###################################################################
        
        """
        Rotate a vector of angle (clockwise)
        
        :param angle: angle
        :param unit: 'dec_deg' or 'radians' (default radians)
        
        :return : rotated components of velocity
        
        """
        
        import numpy as np
        
        if unit == 'degrees':
            angle = np.radians(angle)

        N=self.copy()

        def __rotate__(V,rad_angle):

            cos=np.cos(rad_angle)
            sin=np.sin(rad_angle)
            R=np.array([[cos,-sin],[sin,cos]])
            return(np.dot(R,V))
        

        V=np.array([self.Ve,self.Vn])
        rotated=__rotate__(V,angle)
        (N.Ve,N.Vn)=(rotated[0],rotated[1])
        
        return(N)
 
###################################################################
    def midpoint(self, point , code='Mm'):
###################################################################
        
        """
        returns the mid-point between two GMT_Point
        
        :param point: GMT_Point
        
        :return : mid_point as GMT_Point
        """
        
        import pyacs.lib.coordinates
        
        (Xi,Yi,Zi) = pyacs.lib.coordinates.geo2xyz(self.lon, self.lat, 0.0, 'dec_deg')
        (Xe,Ye,Ze) = pyacs.lib.coordinates.geo2xyz(point.lon, point.lat, 0.0, 'dec_deg')

        X = ( Xi + Xe ) / 2.
        Y = ( Yi + Ye ) / 2.
        Z = ( Zi + Ze ) / 2.
        
        ( lonMm , latMm, heMm ) = pyacs.lib.coordinates.xyz2geo(X, Y, Z, unit='dec_deg' )
        
        
        Mm = GMT_Point(code = code)
        Mm.lon = lonMm
        Mm.lat = latMm        
        Mm.he  = heMm
        
        return Mm

