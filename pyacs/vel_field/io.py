"""I/O methods for Velocity_Field: read and write GMT psvelo files."""


def read(cls, file_name=None, lexclude=[], lonly=[], verbose=False):
    """Read a GMT psvelo file.

    Parameters
    ----------
    file_name : str, optional
        Path to the GMT psvelo file.
    lexclude : list, optional
        List of site codes to exclude.
    lonly : list, optional
        If non-empty, only these site codes are included.
    verbose : bool, optional
        If True, print progress. Default is False.

    Returns
    -------
    Velocity_Field
        Populated velocity field instance.
    """
    import numpy as np

    vf = cls()

    def __gen_fake_code__(n):
        FAKE = []
        for i in np.arange(n):
            fake_code = ("%4s" % hex(i).split('x')[-1].replace('L', '')).replace(' ', '0')
            FAKE.append(fake_code.upper())
        return np.array(FAKE)

    if verbose:
        print("-- Reading GMT psvelo file: %s " % file_name)

    try:
        np_vel = np.atleast_2d(np.genfromtxt(file_name, comments='#'))
    except Exception:
        raise IOError("!!! Could not read file: %s" % file_name)

    if np_vel.size == 0:
        return vf

    if np_vel.shape[1] == 8:
        if verbose:
            print("-- file %s has 8 columns" % file_name)
        np_vel = np.delete(np_vel, -1, axis=1)
        np_code = np.genfromtxt(file_name, comments='#', usecols=(7), dtype=str).reshape(-1)
        code = ''
        np_sort = np.sort(np_code)
        for i in np.arange(np_sort.shape[0]):
            if np_sort[i] == code:
                print("!!!", code)
            code = np_sort[i]
    elif np_vel.shape[1] == 3:
        if verbose:
            print("-- file %s has 3 columns" % file_name)
        np_vel = np.delete(np_vel, -1, axis=1)
        np_code = np.array(np.mat(np.genfromtxt(file_name, comments='#', usecols=(2)))).flatten()
    elif np_vel.shape[1] not in [3, 8]:
        np_code = __gen_fake_code__(np_vel.shape[0])
    else:
        raise IOError("!!! Could not decipher file content: %s", file_name)

    from pyacs.lib.gmtpoint import GMT_Point

    lgmt_points = []
    for i in np.arange(np_vel.shape[0]):
        code = np_code[i]
        if np_vel.shape[1] >= 7:
            lon, lat, Ve, Vn, SVe, SVn, SVen = np_vel[i, :]
            M = GMT_Point(lon=lon, lat=lat, Ve=Ve, Vn=Vn, SVe=SVe, SVn=SVn, SVen=SVen, code=code)
        else:
            lon, lat = np_vel[i, :]
            M = GMT_Point(lon=lon, lat=lat, code=code)

        if verbose:
            M.get_info(display=True)

        if lonly != []:
            if M.code in lonly:
                lgmt_points.append(M)
        else:
            if lexclude != []:
                if M.code not in lexclude:
                    lgmt_points.append(M)
            else:
                lgmt_points.append(M)

    vf.file_name = file_name
    vf.sites = lgmt_points
    return vf


def write(self, out_file=None, lexclude=[], verbose=True, comment='', up=False):
    """Write the velocity field to a GMT psvelo file.

    Parameters
    ----------
    out_file : str, optional
        Output file path. Default is 'tmp_vel.gmt' if None.
    lexclude : list, optional
        List of site codes to exclude from output.
    verbose : bool, optional
        If True, print progress. Default is True.
    comment : str, optional
        Comment line written at top of file.
    up : bool, optional
        If True, write up component instead of east/north. Default is False.
    """
    if out_file is None:
        print("! No file name provided. Using tmp_vel.gmt")
        out_file = 'tmp_vel.gmt'

    if verbose:
        print("-- Writing GMT psvelo file: %s " % out_file)
    fs = open(out_file, 'a')

    if up:
        comment = comment + ' | UP component'

    fs.write('# ' + comment + '\n')

    if up:
        for M in self.sites:
            fs.write("%10.5lf  %10.5lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %s\n" % (M.lon, M.lat, 0.0, M.Vu, 0.0, M.SVu, M.SVen, M.code))
    else:
        for M in self.sites:
            fs.write("%10.5lf  %10.5lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %s\n" % (M.lon, M.lat, M.Ve, M.Vn, M.SVe, M.SVn, M.SVen, M.code))

    fs.close()


def add_point(self, M):
    """Append a GMT_Point to the velocity field.

    Parameters
    ----------
    M : GMT_Point
        Point to add.

    Returns
    -------
    Velocity_Field
        self (for chaining).
    """
    self.sites.append(M)
    return self


def remove_point(self, code):
    """Remove a GMT_Point from the velocity field by code.

    Parameters
    ----------
    code : str
        4-character code of the point to remove.

    Returns
    -------
    Velocity_Field
        New velocity field without the point (self replaced).
    """
    self = self.subset(lexclude=[code])
    return self
