###############################################################################
def geo_to_strike(ilon, ilat, elon, elat):
###############################################################################

    """
    for a given fault segment starting at ilon,lat and ending at elon, elat
    , returns the strike.

    :param ilon,ilat: geographical coordinates of fault segment start point in decimal degrees
    :param elon,elat: geographical coordinates of fault segment end point in decimal degrees
    :returns strike: in decimal degrees clockwise from north

    :note: strike here is taken as the initial bearing. Be cautious with long segments
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
