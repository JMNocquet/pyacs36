"""
Read GipsyX series format time series.
"""


###################################################################
## Loads Gts from GAMIT/GLOBK  pos file format for time series
###################################################################
def read_series(self, tsdir='.', tsfile=None, xyz=True, verbose=False):
    """
    Read GipsyX series time series file.

    Parameters
    ----------
    tsdir : str, optional
        Directory containing the series file.
    tsfile : str, optional
        Series file to load. If None, looks for CODE*.series.
    xyz : bool, optional
        If True, keep/store XYZ; otherwise NEU.
    verbose : bool, optional
        Verbose mode.

    Notes
    -----
    This format does not include the absolute position; relative only.
    """

    # import
    import numpy as np
    import pyacs.lib.astrotime
    import pyacs.lib.coordinates
    import os
    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG


    # name of the file to be read - if not provided, tries to guess

    if (tsfile is None):
        if (self.code is not None):
            from glob import glob
            try:
                pos_file = glob(tsdir + '/' + self.code.upper() + '*.series')[0]
            except:
                ERROR("Could not find any time series file for code ", self.code)
                return ()
        else:
            ERROR("no code or file provided.")

    else:
        pos_file = tsfile

    # file

    self.ifile = os.path.abspath(pos_file)

    # actual read

    # get site code
    from pathlib import Path
    self.code = Path(pos_file).stem
    if len(self.code) != 4:
        WARNING("I was expecting 4-char code for file name. code will be %s" % self.code)

    # reads data
    try:
        self.data = np.atleast_2d( np.genfromtxt(pos_file, usecols=(0,2,1,3,5,4,6,7,9,8)))
    except:
        # there is a big problem
        ERROR("Could not read dN,dE,dU from file: %s " % pos_file, exit=True)

    # check duplicate or non-ordered entries
    if self.data.shape[0] > 1:
        if np.min(np.diff(self.data[:, 0])) <= 0:
            WARNING("time series not properly ordered by dates or dates duplicated. Fixing that.")
            self.reorder()

    # force clean

    self.offsets_dates = []
    self.offsets_values = None
    self.outliers = []
    self.annual = None
    self.semi_annual = None
    self.velocity = None
    self.data_xyz = None

    WARNING("This format is missing longitude, latitude, He information")

    return (self)
