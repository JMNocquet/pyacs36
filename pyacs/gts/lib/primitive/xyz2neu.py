###################################################################
def xyz2neu(self,corr=False,ref_xyz=None, verbose=False):
###################################################################
    """
    populates neu (data) using xyz (data_xyz)
    lon, lat and h will also be set.
    
    :param corr: if True, then standard deviation and correlations will also be calculated  
    :param ref_xyz: [X,Y,Z] corresponding to the 0 of the local NEU frame. If not provided, the first position is used as a reference
    :param verbose: verbose mode
    
    :note: this method is always in place
    """

    import pyacs.lib.coordinates
    import numpy as np

    # Reference point for the local NEU frame
    
    if ref_xyz is not None:
        [xref,yref,zref] = ref_xyz
    else:
        [xref,yref,zref]=self.data_xyz[0,1:4]

    # update X0, Y0, Z0

    self.X0=xref
    self.Y0=yref
    self.Z0=zref
    

    (lam,phi,h) = pyacs.lib.coordinates.xyz2geo(xref,yref,zref)
    R = pyacs.lib.coordinates.mat_rot_general_to_local(lam,phi)

    X0=np.array([xref,yref,zref])
    
    DX=self.data_xyz[:,1:4]-X0
    
    ENU=np.dot(R,DX.T).T
    
    self.data=np.copy(self.data_xyz)
    
    self.data[:,1]=ENU[:,1]
    self.data[:,2]=ENU[:,0]
    self.data[:,3]=ENU[:,2]
    
    self.lon=np.degrees(lam)
    self.lat=np.degrees(phi)
    self.h=h
    
    if corr:
        import pyacs.lib.glinalg
        for i in np.arange(self.data.shape[0]):
            [ _X,_Y,_Z,Sx,Sy,Sz,Rxy,Rxz,Ryz ]=self.data_xyz[i,1:10].tolist()
            CORR_XYZ=np.array([\
                              [1,Rxy,Rxz],\
                              [Rxy,1,Ryz],\
                              [Rxz,Ryz,1]\
                              ])
            STD_XYZ=np.array([Sx,Sy,Sz])
            VCV_XYZ=pyacs.lib.glinalg.corr_to_cov(CORR_XYZ, STD_XYZ)
    
            VCV_ENU=np.dot(np.dot(R,VCV_XYZ),R.T)
            CORR_ENU,STD_ENU=pyacs.lib.glinalg.cov_to_corr(VCV_ENU)
            
            
            self.data[i,4]=STD_ENU[1]
            self.data[i,5]=STD_ENU[0]
            self.data[i,6]=STD_ENU[2]

            self.data[i,7]=CORR_ENU[0,1]
            self.data[i,8]=CORR_ENU[1,2]
            self.data[i,9]=CORR_ENU[0,2]
    
    # reorder and check duplicate dates

    self.reorder(verbose=False)
    # check for duplicates
    #self.correct_duplicated_dates(action='correct',tol= .00000001, in_place=True,verbose=verbose)
 
    return(self)
