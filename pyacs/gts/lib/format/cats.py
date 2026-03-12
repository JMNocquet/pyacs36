"""
Read and write time series files in the CATS format.
"""

###################################################################
def read_cats_file(self,idir='.',ifile=None, gmt=True, verbose=False):
###################################################################
    """
    Read a CATS file and load the time series into this Gts.

    Parameters
    ----------
    idir : str, optional
        Directory for CATS files.
    ifile : str, optional
        Path to CATS file. If None, uses idir/cats_<code>.dat.
    gmt : bool or str, optional
        If True, read lon/lat from ../maps_en_velo.gmt; if str, path to GMT file.
    verbose : bool, optional
        Verbose mode.

    Returns
    -------
    Gts
        self (data loaded from file).
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

    # check that data has 7 columns
    if self.data.shape[1] == 7:
        # adds 3 columns of zeros
        self.data = np.hstack((self.data, np.zeros((self.data.shape[0], 3))))

    return(self)

###################################################################
def write_cats(self,idir='./cats',offsets_dates=None,add_key=''):
###################################################################
    """
    Write a file for CATS analysis.

    If offsets_dates is not None, offset date information is added at the beginning.

    Parameters
    ----------
    idir : str, optional
        Directory to save the file (default './cats').
    offsets_dates : list, optional
        List of offset dates to write in the header.
    add_key : str, optional
        Additional key in the file name (e.g. cats_<code>_<add_key>.dat).

    Returns
    -------
    None
    """

    import numpy as np
    import os

    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG


    if not os.path.isdir(idir):
        try:
            os.makedirs( idir, exist_ok=True)
        except:
            WARNING("Could not create directory ",idir)

    
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
