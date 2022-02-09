###################################################################
## Loads Gts from tdp (Gipsy kinmeatics provided by Cedric Twardzik 17/04/2018) format for time series
###################################################################
def read_tdp(self,idir='.',ifile=None, gmt=True, verbose=False):
    """
    Read tdp (Gipsy kinematics provided by Cedric Twardzik 17/04/2018) format for time series
    """
    
    import numpy as np
    import os
    import pyacs.lib.astrotime
    
    # file
    self.ifile=os.path.abspath(ifile)
        
    raw_data=np.genfromtxt(ifile)
    
    # converts dates to decimal year
    
    np_decyear = pyacs.lib.astrotime.datetime2decyear(pyacs.lib.astrotime.datetime_from_calarray(raw_data[:,:6]))
    
    # creates empty data attribute
    self.data=np.zeros((raw_data.shape[0],10)) + 1.E-3
    
    # populates data
    self.data[:,0] = np_decyear
    self.data[:,1:4]=raw_data[:,7:10]
