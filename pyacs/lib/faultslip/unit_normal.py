def unit_normal( strike , dip ):
    """
    returns the unit vector normal to the fault plane in ENU (Up positive upward)

    :param strike: fault strike in degrees
    :param dip: fault dip in degrees
    :return:
    """

    import numpy as np

    rstrike = np.radians( strike )
    rdip    = np.radians( dip )

    unit_E =  np.sin( rdip ) * np.cos( rstrike )

    unit_N = -np.sin( rdip ) * np.sin( rstrike )

    unit_U = np.cos( rdip )

    return( unit_E  , unit_N , unit_U )