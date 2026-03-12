###################################################################
## change reference frame for the time series
###################################################################

def frame(self, frame=None, verbose=False):
    """
    Rotates a time series according to an Euler pole.
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
        ERROR(error)
        return (self)
    ###########################################################################

    # Euler poles taken from pygvel_pole_info.py
    lEuler = {}
    lEuler['soam'] = [-132.21, -18.83, 0.121]
    lEuler['nas'] = [-97.52, 6.48, 0.359]
    lEuler['nazca'] = [-94.4, 61.0, 0.57]
    lEuler['nas_wrt_soam'] = [-83.40, 15.21, 0.287]
    lEuler['inca_wrt_soam'] = [-63.76, 22.47, 0.092]
    lEuler['eura'] = [-98.83483039304333, 54.22539546556553,
                      0.25678223107826376]  # from Altamimi et al., 2012, eura_wrt_itrf2008

    euler_vector = np.array(lEuler[frame])

    new_Gts = self.remove_pole(euler_vector, verbose=verbose)
    return (new_Gts)
