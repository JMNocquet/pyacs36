###################################################################
## remove_pole
###################################################################

def remove_pole(self, pole, pole_type='euler', verbose=True):
    """
    Remove velocity predicted by an Euler pole or a rotation rate vector from a time series.
    pole is a 1D array with 3 values.
    Requires self.lon & self.lat attributes to have been filled before.
    Returns a new Gts instance.
    """

    import numpy as np
    from pyacs.gts.Gts import Gts
    import inspect
    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG

    # after this method .data  and .data_xyz are not consistent so .data_xyz is set to None
    #self.data_xyz = None


    ###########################################################################
    # check data is not None
    from pyacs.gts.lib.errors import GtsInputDataNone

    try:
        if self.data is None:
            # raise exception
            raise GtsInputDataNone(inspect.stack()[0][3], __name__, self)
    except GtsInputDataNone as error:
        # print PYACS WARNING
        WARNING(error)
        return (self)
    ###########################################################################

    if (self.lon is None) or (self.lat is None):
        ERROR("lon & lat needs to be properly filled to use this method.")
        return ()

    from pyacs.lib.gmtpoint import GMT_Point

    M = GMT_Point(code=self.code, lon=self.lon, lat=self.lat, Ve=0., Vn=0., SVe=0., SVn=0.)
    # N=M.substract_pole(pole,pole_type)

    N = M.pole(W=pole, SW=None, type_euler='euler', option='predict')

    vel_neu = np.array([N.Vn, N.Ve, 0.]) * 1E-3

    VERBOSE("Removing velocity (NEU): %10.3lf %10.3lf %10.3lf" % tuple((vel_neu * 1.E3).tolist()))

    new_Gts = self.remove_velocity(vel_neu)
    return (new_Gts)
