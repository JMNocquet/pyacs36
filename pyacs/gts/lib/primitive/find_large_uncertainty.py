def find_large_uncertainty( self , sigma_thresold=10 , verbose=True, lcomponent='NE' ):
    """
    Find dates with large uncertainty and flag them as outliers.
    
    :param sigma_threshold: value (mm) for a date to be flagged.
    :param verbose: verbose mode
    :param lcomponent: list of components to be checked. default = 'NE'
    """
    
    # import
    import numpy as np

    # lindex
    lidxn = []
    lidxe = []
    lidxu = []
    
    # threshold in mm
    stmm = sigma_thresold * 1.E-3
    
    # North component
    if 'N' in lcomponent:
        lidxn = np.argwhere( self.data[:,4] > stmm ).flatten().tolist()
    
    # East component
    if 'E' in lcomponent:
        lidxe = np.argwhere( self.data[:,5] > stmm ).flatten().tolist() 

    # Up component
    if 'U' in lcomponent:
        lidxu = np.argwhere( self.data[:,6] > stmm ).flatten().tolist()
    
    # update
    
    self.outliers = list(set(sorted(self.outliers + lidxn + lidxe + lidxu)))
    
    return self