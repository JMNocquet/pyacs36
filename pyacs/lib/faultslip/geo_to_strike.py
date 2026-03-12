###############################################################################
def geo_to_strike(ilon, ilat, elon, elat):
###############################################################################

    """Return strike for a fault segment from start and end coordinates.

    Parameters
    ----------
    ilon : float
        Longitude of fault segment start in decimal degrees.
    ilat : float
        Latitude of fault segment start in decimal degrees.
    elon : float
        Longitude of fault segment end in decimal degrees.
    elat : float
        Latitude of fault segment end in decimal degrees.

    Returns
    -------
    float
        Strike in decimal degrees clockwise from north.

    Notes
    -----
    Strike is the initial bearing; be cautious with long segments.
    """

    import numpy as np

    # radian conversion
    [lon1, lat1, lon2, lat2] = np.radians([ilon, ilat, elon, elat])

    strike_rad_from_east_ccwise = np.arctan2(
        np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(lon2 - lon1), \
        np.sin(lon2 - lon1) * np.cos(lat2))

    strike_deg_from_north_cwise = 90. - np.degrees(strike_rad_from_east_ccwise)

    if strike_deg_from_north_cwise > 180.:
        strike_deg_from_north_cwise = strike_deg_from_north_cwise - 360.

    return (strike_deg_from_north_cwise)
