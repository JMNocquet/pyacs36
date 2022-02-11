###################################################################
## change reference frame for the time series
###################################################################

def frame(self, frame=None, in_place=False, verbose=False):
    """
    Rotates a time series according to an Euler pole
    Returns a new Gts instance
    """

    import numpy as np
    from pyacs.gts.Gts import Gts
    import inspect

    # after this method .data  and .data_xyz are not consistent so .data_xyz is set to None
    self.data_xyz = None

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
    if in_place:
        self.data = new_Gts.data.copy()
    return (new_Gts)
