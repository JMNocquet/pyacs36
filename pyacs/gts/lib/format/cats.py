"""
Reads and writes time series files in the CATS format
"""

###################################################################
def read_cats_file(self,idir='.',ifile=None, gmt=True, verbose=False):
###################################################################
    """
    Read cats files in a directory and actually loads the time series
    
    """
    
    import numpy as np
    import os
    
    if gmt==True:
        gmt_file=idir+'/../maps_en_velo.gmt'
    if  isinstance(gmt,str):
        gmt_file=gmt
    
    if gmt != False:
        self.read_lon_lat(gmt_file,verbose=verbose)
    
    if ifile is None:
        cats_file_basename= idir + '/cats_'+self.code+'.dat'
    else:
        cats_file_basename=ifile

    # file
    self.ifile=os.path.abspath(cats_file_basename)
    
    self.data=np.genfromtxt(cats_file_basename,comments='#')

###################################################################
def write_cats(self,idir,offsets_dates=None,add_key=''):
###################################################################
    """
    Writes a file for a cats analysis
    if offsets_dates is not None then offsets are added at the beginning of the file
    """

    import numpy as np
    import os

    if not os.path.isdir(idir):
        try:
            os.mkdir(idir)
            print("-- Creating directory ",idir)
        except:
            print("!!! Error : Could not create directory ",idir)
    

    
    if add_key != '':
        cats_file=idir+'/cats_'+self.code+'_'+add_key+'.dat'
    else:
        cats_file=idir+'/cats_'+self.code+'.dat'
    
    # header lines for outliers
    header=''
    if self.offsets_dates is not None:
        for date in self.offsets_dates:header=header+("%10.5lf 7\n" % date)

    if offsets_dates is not None:
        for date in offsets_dates:header=header+("%10.5lf 7\n" % date)
            
    
    # main data
    np.savetxt(cats_file, self.data[:,:7], fmt="%10.5lf %10.5lf %10.5lf %10.5lf %10.5lf %10.5lf %10.5lf",header=header)

    return(self)
