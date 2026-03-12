def unit_normal( strike , dip ):
    """Return unit vector normal to the fault plane in ENU (Up positive upward).

    Parameters
    ----------
    strike : float
        Fault strike in degrees.
    dip : float
        Fault dip in degrees.

    Returns
    -------
    tuple of float
        (unit_E, unit_N, unit_U) in ENU.
    """

    import numpy as np

    rstrike = np.radians( strike )
    rdip    = np.radians( dip )

    unit_E =  np.sin( rdip ) * np.cos( rstrike )

    unit_N = -np.sin( rdip ) * np.sin( rstrike )

    unit_U = np.cos( rdip )

    return( unit_E  , unit_N , unit_U )