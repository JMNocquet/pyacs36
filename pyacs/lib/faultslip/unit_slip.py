def unit_slip(strike, dip, rake):
    """Return unit slip vector in ENU convention (Up positive).

    Parameters
    ----------
    strike, dip, rake : float
        Fault angles in decimal degrees.

    Returns
    -------
    tuple
        (unit_E, unit_N, unit_U) unit slip components in ENU.
    """

    import numpy as np

    rstrike = np.radians( strike )
    rrake   = np.radians( rake )
    rdip    = np.radians( dip )

    unit_E = np.cos( rrake )* np.sin( rstrike ) - np.cos( rdip ) * np.sin( rrake ) * np.cos( rstrike )

    unit_N = np.cos( rrake )* np.cos( rstrike ) + np.cos( rdip ) * np.sin( rrake ) * np.sin( rstrike )

    unit_U = - np.sin( rdip ) * np.sin( rrake )

    return( unit_E, unit_N, unit_U )