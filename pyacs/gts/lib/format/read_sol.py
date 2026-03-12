"""Read SGC sol format for time series."""


###################################################################
## Loads Gts from GAMIT/GLOBK  pos file format for time series
###################################################################
def read_sol(self,tsdir='.',ifile=None, verbose=False):
    """Read SGC sol file and load the time series into this Gts instance.

    Parameters
    ----------
    tsdir : str, optional
        Directory of sol/pos files. Default is '.'.
    ifile : str, optional
        Path to the sol file to load. If None, a file CODE*.pos may be used.
    verbose : bool, optional
        If True, print progress. Default is False.

    Returns
    -------
    Gts
        self (data, code, X0, Y0, Z0, t0 are populated from file).

    Notes
    -----
    A pos/sol file contains (almost) all needed info. If ifile is None,
    read_pos looks for a file named CODE*.pos.
    """

    # import
    import numpy as np
    import pyacs.lib.astrotime as at
    import pyacs.lib.coordinates
    import os
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG
    from datetime import datetime
    from pyacs.gts.Gts import Gts
    
    # check whether the file exists
    if not os.path.exists(ifile):
        ERROR("File %s does not exist" % ifile, exit=True)

    # read the file
    # read date from first column
    dates=np.genfromtxt(ifile,usecols=(0),dtype=str)
    # convert dates from '2013-06-14' to decyear using datetime.strptime
    dates_decyear = at.datetime2decyear([datetime.strptime(date, '%Y-%m-%d') for date in dates])
    # reads xyz_sxyz_corr_xyz
    data_xyz=np.genfromtxt(ifile,usecols=(0,2,3,4,5,6,7,8,9,10),dtype=float)
    # data[4:8] are the standard deviations and data[8:11] are the correlations
    # add dates_decyear to data_xyz
    data_xyz[:,0]=dates_decyear
    # create a Gts object
    self.ifile=ifile
    self.code=os.path.basename(ifile).split('.')[0]
    self.data_xyz = data_xyz
    self.xyz2neu( corr=True)

    return self