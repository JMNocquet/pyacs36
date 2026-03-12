###################################################################
def sel_from_grid(self, grid, depth_range=[0,200]):
###################################################################
    """Select sites above a grid (e.g. slab2) within a depth range.

    Parameters
    ----------
    grid : str
        Path to NetCDF grid file.
    depth_range : list, optional
        [depth_min_km, depth_max_km]. Default is [0, 200].

    Returns
    -------
    Sgts
        New Sgts with selected time series.
    """

    # import
    from pyacs.gts.Sgts import Sgts

    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG
    import pyacs.debug

    import numpy as np
    from tqdm import tqdm
    import scipy.interpolate
    import netCDF4



    import inspect

    VERBOSE("Running Sgts.%s" % inspect.currentframe().f_code.co_name)

    # reads the grid

    try:
        ds = netCDF4.Dataset(grid)
    except:
        ERROR(("Could not read grid: %s" % (grid)), exit=True)

    # case grid -180/180 or 0/360
    longitudes = ds.variables['x'][:]
    if np.max(longitudes) > 180.:
        case_grid_type = '0_360'
    else:
        case_grid_type = '-180_180'

    # creates the interpolator
    z = ds.variables['z'][:].T
    # change JMN 28/03/2023 - account for 0-360 grid like slab2
    #longitudes = ds.variables['x'][:]
    #lidx = np.where( longitudes > 180 )[0]
    #if lidx.shape[0] > 0:
    #    longitudes = longitudes - 360.
    #longitudes[lidx] = longitudes[lidx] - 360.
    interp = scipy.interpolate.RegularGridInterpolator(
        tuple((longitudes, ds.variables['y'][:])),
        z.data,
        method='linear',
        bounds_error=False)

#    interp = scipy.interpolate.RegularGridInterpolator(
#        tuple((ds.variables['x'][:], ds.variables['y'][:])),
#        z.data,
#        method='linear',
#        bounds_error=False)

    # initialize new Sgts
    new_Sgts = Sgts(read=False)

    # create array of lon,lat,depth
    COOR = np.zeros(( self.n() , 3 ))
    np_name = np.array(self.lcode(),dtype=str)
    for i in np.arange( np_name.shape[0] ):
        COOR[i,0] = self.__dict__[ np_name[i] ].lon
        COOR[i,1] = self.__dict__[ np_name[i] ].lat

    # case_grid_type == '0_360':
    if case_grid_type == '0_360':
        lidx = np.where( COOR[:,0] < 0 )[0]
        COOR[lidx,0] = COOR[lidx,0] + 360.
    if case_grid_type == '-180_180':
        lidx = np.where( COOR[:,0] > 180 )[0]
        COOR[lidx,0] = COOR[lidx,0] - 360.

    # get depth
    COOR[:,2] = np.fabs( interp( (COOR[:,0], COOR[:,1]) ) )

    # get index
    lidx = np.where( (COOR[:,2] > depth_range[0]) & (COOR[:,2] < depth_range[1]) )[0]
    VERBOSE("#sites selected: %d / %d" % (len(lidx),np_name.shape[0]) )

    #return
    return self.sub(linclude=np_name[lidx])
