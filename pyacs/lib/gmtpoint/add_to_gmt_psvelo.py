"""Add GMT_Point to a GMT psvelo file."""


def add_to_gmt_psvelo(self, fname, overwrite=False, verbose=False):
    """Add the current GMT_Point to a GMT psvelo file.

    Parameters
    ----------
    fname : str
        Path to the GMT psvelo file (created if missing).
    overwrite : bool, optional
        If True, replace existing line with same code. Default is False.
    verbose : bool, optional
        If True, print progress. Default is False.
    """
    import pyacs.vel_field
    vel = pyacs.vel_field.Velocity_Field()
    try:
        vel.read(file_name=fname, verbose=verbose)
    except Exception:
        if verbose:
            print('-- file does not exist. Creating: %s ', fname)

    if overwrite:
        if self.code in vel.lcode():
            if verbose:
                print('-- replacing site %s:%s' % (self.code, vel.site(self.code).get_info()))
                print('-- with      site %s:%s' % (self.code, self.get_info()))

        vel.remove_point(self.code)

    vel.add_point(self)

    vel.write(fname, verbose=verbose)
