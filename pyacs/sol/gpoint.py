class Gpoint:
    
    """
    Geodetic point
    
    Used with sinex class
    
    """

    def __init__ (self,X=None,Y=None,Z=None,\
                  SX=None,SY=None,SZ=None,
                  epoch=None,code=None,pt=None,soln=None,domes=None,\
                  VX=None,VY=None,VZ=None,\
                  SVX=None,SVY=None,SVZ=None,\
                  Ve=None,Vn=None,Vu=None,\
                  SVe=None,SVn=None,SVu=None,\
                  lon=None,lat=None,he=None,\
                  Cv_xyz=None, Cv_enu=None, index=None):
        self.X=X
        self.Y=Y
        self.Z=Z
        self.SX=SX
        self.SY=SY
        self.SZ=SZ
        self.epoch=epoch
        self.code=code
        self.pt=pt
        self.soln=soln
        self.domes=domes
        self.VX=VX
        self.VY=VY
        self.VZ=VZ
        self.SVX=SVX
        self.SVY=SVY
        self.SVZ=SVZ
        self.Ve=Ve
        self.Vn=Vn
        self.Vu=Vu
        self.SVe=SVe
        self.SVn=SVn
        self.SVu=SVu
        self.lon=lon
        self.lat=lat
        self.he=he
        self.Cv_xyz=Cv_xyz
        self.Cv_enu=Cv_enu
        self.index=index
        
    def assign_index(self,index):
        self.index=index
 #   def get_index(self):
 #       return self.index
        
        
    def posxyz(self):
        return(self.X,self.Y,self.Z)
    
    def staxyz(self):
        return (self.SX, self.SY, self.SZ)
    
    def covarxyz(self):
        return self.Cv_xyz

    def covarenu(self):
        return self.Cv_enu

    def velxyz(self):
        return(self.VX,self.VY,self.VZ)

    def read_from_sinex(self,sinex_name):
        """Reads an populate self from a sinex file"""
        import pyacs.sol.sinex as SSinex
        sinex=SSinex.SSinex(sinex_name)
        sinex.read()
        M=sinex.subset(lcode=[self.code])[0]
        self.X=M.X
        self.Y=M.Y
        self.Z=M.Z
        self.SX=M.SX
        self.SY=M.SY
        self.SZ=M.SZ

        self.VX=M.VX
        self.VY=M.VY
        self.VZ=M.VZ
        self.SVX=M.SVX
        self.SVY=M.SVY
        self.SVZ=M.SVZ

        self.epoch=M.epoch
        self.pt=M.pt
        self.soln=M.soln
        self.domes=M.domes
        
    def write_as_estimate(self,fname,index=1,VEL=True, VENU=None,CVENU=None,ENU=None,soln=1,epoch=2010.0):
        """Writes site information formatted for a ESTIMATE section of a sinex file"""
#   139 STAX   GENO  A    1 10:029:43185 m    2 0.450789227114383E+07 .204988E+00
        import pyacs.lib.astrotime as AstroTime
        import pyacs.lib.coordinates as Coordinates
        import numpy as np
        if VENU is not None:
            pass
        
        def decyear2yeardoysod(decyear):
            """ 
                year,doy,sod = decyear2yeardoysod(decyear)
                Converts decimal year into calendar year, day of year, second of day
            """
            year=int(decyear)
            if year >= 2000: year = year - 2000
            else: year = year - 1900
            
            doy,ut = AstroTime.decyear2dayno(decyear)
            
            sod = int(ut*86400)
            if sod !=0:
                doy += 1
                sod = 0
            return year,doy,sod

        year,doy,sod=decyear2yeardoysod(epoch)


        (lam,phi,h)=Coordinates.xyz2geo(self.X,self.Y,self.Z)
        R=Coordinates.mat_rot_local_to_general(lam,phi,h)
        
        if ENU  is not None:
            delta_XYZ=np.dot(R,ENU.reshape(3,1))
            self.X=self.X-delta_XYZ[0]/1000.0
            self.Y=self.Y-delta_XYZ[1]/1000.0
            self.Z=self.Z-delta_XYZ[2]/1000.0
        
        if VENU!=None:
            VXYZ=np.dot(R,VENU.reshape(3,1))
            #print VXYZ*1000.0
             
            self.VX,self.VY,self.VZ=VXYZ[0,0],VXYZ[1,0],VXYZ[2,0]
            (self.SVX,self.SVY,self.SVZ)=(1.0,1.0,1.0)
            self.X=self.X+(self.epoch-epoch)*self.VX
            self.Y=self.Y+(self.epoch-epoch)*self.VY
            self.Z=self.Z+(self.epoch-epoch)*self.VZ
        else:
            if self.VX==None:
                (self.VX,self.VY,self.VZ,self.SVX,self.SVY,self.SVZ)=(0.0,0.0,0.0,1.0,1.0,1.0)
            if CVENU!=None:
                CVXYZ=np.dot(R,CVENU.reshape(3,1))
                #print VXYZ*1000.0
                 
                self.VX,self.VY,self.VZ=self.VX-CVXYZ[0,0],self.VY-CVXYZ[1,0],self.VZ-CVXYZ[2,0]
                (self.SVX,self.SVY,self.SVZ)=(1.0,1.0,1.0)
                self.X=self.X+(self.epoch-epoch)*self.VX
                self.Y=self.Y+(self.epoch-epoch)*self.VY
                self.Z=self.Z+(self.epoch-epoch)*self.VZ
                
        
        #import sys
        #sys.exit()
        
        fs=open(fname,'w')
        print('index ',index)
        fs.write("%6d STAX   %4s  A %04d %02d:%03d:%05d m    2 %21.14E %11.5E\n" % (index,self.code,soln,year,doy,sod,self.X,self.SX))
        index=index+1
        fs.write("%6d STAY   %4s  A %04d %02d:%03d:%05d m    2 %21.14E %11.5E\n" % (index,self.code,soln,year,doy,sod,self.Y,self.SY))
        index=index+1
        fs.write("%6d STAZ   %4s  A %04d %02d:%03d:%05d m    2 %21.14E %11.5E\n" % (index,self.code,soln,year,doy,sod,self.Z,self.SZ))
        if not VEL:
            fs.close()
            return()
        index=index+1
        fs.write("%6d VELX   %4s  A %04d %02d:%03d:%05d m    2 %21.14E %11.5E\n" % (index,self.code,soln,year,doy,sod,self.VX,self.SVX))
        index=index+1
        fs.write("%6d VELY   %4s  A %04d %02d:%03d:%05d m    2 %21.14E %11.5E\n" % (index,self.code,soln,year,doy,sod,self.VY,self.SVY))
        index=index+1
        fs.write("%6d VELZ   %4s  A %04d %02d:%03d:%05d m    2 %21.14E %11.5E\n" % (index,self.code,soln,year,doy,sod,self.VZ,self.SVZ))

        fs.close()
        
        return()
    
    def apply_helmert(self,T):
        
        import pyacs.sol.helmert as Helmert
        import copy
        from numpy import array
        
        H=Helmert.helmert_matrix(self.X,self.Y,self.Z,tran=True,rot=True,scale=True,equilibrate=True)
        
        X1=array([self.posxyz()])

        X=H*T+X1.T
        M=copy.deepcopy(self)
        M.X=X[0,0]
        M.Y=X[1,0]
        M.Z=X[2,0]
        return M


    def substract_pole(self,W=None,type='rot'):
        """
        Substract the velocity predicted by an Euler pole or a rotation rate vector
        """
        import numpy as np
        if isinstance(W,list):W=np.array(W)
        
        if len(W.shape)==1:W=W.reshape(3,1)


        import pyacs.lib.coordinates as Coordinates
        import pyacs.lib.euler as euler
        import numpy as np

        ltype=['euler','rot'] # type is either Euler pole (long. lat in dec. degrees and deg/Myr or cartesian rotation rate vector in rad/yr)
        if (type not in ltype):
            print("==> type should be either 'euler' or 'rot' ")
        if (type == 'euler'):
            # convert to rotation rate vector
            (lon, lat, omega)=(W[0,0],W[1,0],W[2,0])
            (wx,wy,wz)=euler.euler2rot(lon, lat, omega)
            W=np.zeros([3,1])
            W[0,0]=wx
            W[1,0]=wy
            W[2,0]=wz

        import math

        

        R=Coordinates.mat_rot_general_to_local(math.radians(self.lon),math.radians(self.lat),0.)
        (x,y,z)=Coordinates.geo2xyz(math.radians(self.lon),math.radians(self.lat),0.)
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
        
        import copy
        N=copy.deepcopy(self)
        
        N.Ve=self.Ve-pve
        N.Vn=self.Vn-pvn

        return(N)
        




    def xyz_distance(self,M):
        from math import sqrt
        return(sqrt((self.X-M.X)**2+(self.Y-M.Y)**2+(self.Z-M.Z)**2))
