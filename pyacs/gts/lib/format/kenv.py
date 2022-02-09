"""
Reads time series files in the kenv format
"""

###################################################################
## Loads Gts from kenv (magnet) format for time series
###################################################################
def read_kenv(self,ifile,date_type='jd'):
    """
    Read kenv file (magnet) format for time series
    """
    
    import numpy as np
    import os

    
    # file
    self.ifile=os.path.abspath(ifile)
        
    raw_data=np.genfromtxt(ifile,skip_header=1)
    self.data=np.zeros((raw_data.shape[0],10))
    # deal with dates
    if date_type=='jd':
        self.data[:,0]=raw_data[:,6]+raw_data[:,7]/86400
    self.data[:,1:4]=raw_data[:,8:11]
    self.data[:,5:]=raw_data[:,14:18]

