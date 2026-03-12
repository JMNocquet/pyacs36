###################################################################
def save_velocity(self,vel_file='../stat/vel', en=True, up=True):
###################################################################
    """Save horizontal and up velocities in GMT psvelo format from .velocity attributes.

    Parameters
    ----------
    vel_file : str, optional
        Output basename (e.g. '../stat/vel' -> _en.gmt, _up.gmt). Default is '../stat/vel'.
    en : bool, optional
        If True, write East & North velocity. Default is True.
    up : bool, optional
        If True, write Up velocity. Default is True.

    Returns
    -------
    Sgts
        self.
    """

    # import
    import numpy as np
    import pyacs.lib.utils


    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG

    import inspect

    VERBOSE("Running Sgts.%s" % inspect.currentframe().f_code.co_name)

    # output file name
    h_vel=vel_file+'_en.gmt'
    u_vel=vel_file+'_up.gmt'

    # initialize
    np_name = np.array(self.lcode() , dtype=str )
    np_vel_en = np.zeros(( np_name.shape[0] , 7 ))
    np_vel_en[:,3] =np.nan
    np_vel_up = np.zeros((np_name.shape[0], 7))
    np_vel_up[:, 3] = np.nan

    # loop on sites (en)
    for i in np.arange( np_name.shape[0] ):
        code = np_name[i]
        np_vel_en[i,0] = self.__dict__[code].lon
        np_vel_en[i,1] = self.__dict__[code].lat
        np_vel_up[i,0] = self.__dict__[code].lon
        np_vel_up[i,1] = self.__dict__[code].lat
        try:
            np_vel_en[i, 2] = self.__dict__[code].velocity[1]*1.E3
            np_vel_en[i, 3] = self.__dict__[code].velocity[0]*1.E3
            np_vel_en[i, 4] = self.__dict__[code].velocity[4]*1.E3
            np_vel_en[i, 5] = self.__dict__[code].velocity[3]*1.E3

            np_vel_up[i, 3] = self.__dict__[code].velocity[2]*1.E3
            np_vel_up[i, 5] = self.__dict__[code].velocity[5]*1.E3

        except:
            WARNING( ("No velocity for %s" % code) )

    # remove nan
    np_index_nan = np.where( np.isnan(np_vel_en[:,3]) )[0]
    if np_index_nan.shape[0] != 0:
        WARNING(("Removing %d sites from output" % np_index_nan.shape[0]))
        np_vel_en = np.delete(np_vel_en,np_index_nan,axis=0)
        np_name = np.delete(np_vel_en,np_index_nan)

    # save
    fmt = "%10.5lf  %10.5lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %s"
    if en:
        MESSAGE(("saving %s with %d sites" % (h_vel,np_name.shape[0]) ))
        pyacs.lib.utils.save_np_array_with_string(np_vel_en,np_name,fmt,h_vel,comment="# lon lat ve vn sve svn sven code - mm/yr")

    if up:
        MESSAGE(("saving %s with %d sites" % (u_vel, np_name.shape[0])))
        pyacs.lib.utils.save_np_array_with_string(np_vel_up, np_name, fmt, u_vel,
                                                  comment="# lon lat 0 vu 0 svu 0 code - mm/yr")

    return(self)
