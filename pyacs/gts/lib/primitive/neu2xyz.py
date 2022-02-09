###################################################################
def neu2xyz(self,corr=False,verbose=False):
###################################################################
    """
    populates .data_xyz from .data
    requires X0,Y0,Z0 attributes to be set
    
    :param corr: if True, then standard deviation and correlations will also be calculated  
    :param verbose: verbose mode
    
    """

    # check X0,Y0,Z0
    if self.X0 is None:
        print("!!! ERROR: X0,Y0,Z0 attributes required for neu2xyz method.")
        self.data_xyz = None
        return(self)
    
    # import

    import pyacs.lib.coordinates
    import numpy as np

    # reference pos
    
    xref = self.X0
    yref = self.Y0
    zref = self.Z0
    #X0=np.array([xref,yref,zref])

    # local frame to geocentric frame rotation matrix - assumes ENU convention
     
    (lam,phi,_h) = pyacs.lib.coordinates.xyz2geo(xref,yref,zref)
    R = pyacs.lib.coordinates.mat_rot_local_to_general(lam, phi)

    DNEU =  np.copy(self.data[:,1:4])
    # swap for ENU convention - pyacs.gts is NEU
    DENU =  DNEU
    DENU[:,[0, 1]] = DENU[:,[1, 0]]
    
    # DXYZ
    
    DXYZ = np.dot(R,DENU.T).T
    
    self.data_xyz = np.zeros( ( self.data.shape[0] , 10 ) )
    self.data_xyz[:,0] = self.data[:,0].copy()
    
    self.data_xyz[:,1]=DXYZ[:,0] + self.X0
    self.data_xyz[:,2]=DXYZ[:,1] + self.Y0
    self.data_xyz[:,3]=DXYZ[:,2] + self.Z0
    
    
    if corr:
        import pyacs.lib.glinalg
        for i in np.arange(self.data.shape[0]):
            [  _dN , _dE , _dU , Sn , Se , Su , Rne , Rnu , Reu ]=self.data[i,1:10].tolist()
            CORR_ENU=np.array([\
                              [1,Rne,Reu],\
                              [Rne,1,Rnu],\
                              [Reu,Rnu,1]\
                              ])
            STD_ENU=np.array([Se,Sn,Su])
            VCV_ENU=pyacs.lib.glinalg.corr_to_cov(CORR_ENU, STD_ENU)
    
            VCV_XYZ=np.dot(np.dot(R,VCV_ENU),R.T)
            CORR_XYZ,STD_XYZ=pyacs.lib.glinalg.cov_to_corr(VCV_XYZ)
            
            self.data_xyz[i,4]=STD_XYZ[0]
            self.data_xyz[i,5]=STD_XYZ[1]
            self.data_xyz[i,6]=STD_XYZ[2]

            self.data_xyz[i,7]=CORR_XYZ[0,1]
            self.data_xyz[i,8]=CORR_XYZ[0,2]
            self.data_xyz[i,9]=CORR_XYZ[1,2]
    
    # reorder and check duplicate dates

    self.reorder(verbose=False)
    # check for duplicates
    self.correct_duplicated_dates(action='correct',tol= .00000001, in_place=True,verbose=verbose)
 
    return(self)
