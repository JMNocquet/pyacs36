"""Print information for GMT_Point."""


def get_info(self, display=False, legend=False, verbose=False):
    """Print information for the current GMT_Point.

    Parameters
    ----------
    display : bool, optional
        If True, print to stdout. Default is False.
    legend : bool, optional
        If True, print column legend. Default is False.
    verbose : bool, optional
        If True, print verbose message. Default is False.

    Returns
    -------
    str
        Formatted info string.
    """
    info_legend = "code     "
    info = ("%4s " % self.code)

    info_legend = info_legend + " long.        lat."
    info_legend = info_legend + "       Ve.        Vn."
    info_legend = info_legend + "      SVe.       SVn.       SVen"
    info_legend = info_legend + "      V_magnitud    V_azimuth"

    info = info + ("%10.5lf  %10.5lf" % (self.lon, self.lat))

    if self.Ve is not None:
        info = info + ("%10.3lf %10.3lf" % (self.Ve, self.Vn))
        (mag, az) = self.magaz()

        if self.SVe is not None:
            info = info + ("%10.3lf %10.3lf %10.3lf" % (self.SVe, self.SVn, self.SVen))
            info = info + ("      %10.3lf   %10.1lf " % (mag, az))

    if legend:
        print(info_legend)
        print(info)

    if display:
        print(info)

    if verbose:
        print('-- returning info string for site: %s', self.code)

    return info
