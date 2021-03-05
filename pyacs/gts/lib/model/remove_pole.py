###################################################################
## remove_pole
###################################################################

def remove_pole(self, pole, pole_type='euler', in_place=False, verbose=True):
    """
    remove velocity predicted by an Euler pole or a rotation rate vector from a time series
    pole is a 1D array with 3 values
    requires self.lon & self.lat attributes to have been filled before
    if in_place = True then replace the current time series
    """

    import numpy as np
    from pyacs.gts.Gts import Gts
    import inspect


    ###########################################################################
    # check data is not None
    from pyacs.gts.lib.errors import GtsInputDataNone

    try:
        if self.data is None:
            # raise exception
            raise GtsInputDataNone(inspect.stack()[0][3], __name__, self)
    except GtsInputDataNone as error:
        # print PYACS WARNING
        print(error)
        return (self)
    ###########################################################################

    if (self.lon is None) or (self.lat is None):
        print("!!! ERROR: lon & lat needs to be properly filled to use this method.")
        return ()

    from pyacs.lib.gmtpoint import GMT_Point

    M = GMT_Point(code=self.code, lon=self.lon, lat=self.lat, Ve=0., Vn=0., SVe=0., SVn=0.)
    # N=M.substract_pole(pole,pole_type)

    N = M.pole(W=pole, SW=None, type_euler='euler', option='predict')

    vel_neu = np.array([N.Vn, N.Ve, 0.]) * 1E-3

    if verbose: print("-- Removing velocity (NEU)", vel_neu * 1.E3)

    new_Gts = self.remove_velocity(vel_neu)

    if in_place:
        self.data = new_Gts.data.copy()
    return (new_Gts)
