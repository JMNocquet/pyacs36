"""Info and accessor methods for Velocity_Field."""

import numpy as np


def info(self, details=True):
    """Print basic information about the velocity field.

    Parameters
    ----------
    details : bool, optional
        If True, list all site codes. Default is True.
    """
    print("-- velocity field from: %s " % self.file_name)
    print("-- number of sites    : %d " % len(self.sites))

    if details:
        n_site_per_line = 10
        print('-- list of sites: ')
        A = np.append(np.array(self.lcode()), ['    '] * (n_site_per_line - np.remainder(len(self.lcode()), n_site_per_line))).reshape(-1, n_site_per_line)
        print(np.array2string(A, separator='').replace('[', ' ').replace(']', ' ').replace('\'', ' '))


def nsites(self):
    """Return the number of sites in the velocity field.

    Returns
    -------
    int
        Number of sites.
    """
    return len(self.sites)


def l_GMT_Point(self):
    """Return the velocity field as a list of GMT_Point objects.

    Returns
    -------
    list
        List of GMT_Point instances.
    """
    lsite = []
    for M in self.sites:
        lsite.append(M)
    return lsite


def print_info_site(self, code, verbose=False):
    """Print information for a site given its code.

    Parameters
    ----------
    code : str
        4-character site code.
    verbose : bool, optional
        If True, print legend. Default is False.
    """
    for M in self.sites:
        if M.code == code:
            M.get_info(legend=verbose, display=True)


def lcode(self):
    """Return a list of all point codes in the velocity field.

    Returns
    -------
    list
        List of 4-character site codes.
    """
    lcode_out = []
    for M in self.sites:
        if M.code not in lcode_out:
            lcode_out.append(M.code)
    return lcode_out


def site(self, code):
    """Return a site as a GMT_Point by code.

    Parameters
    ----------
    code : str
        4-character site code.

    Returns
    -------
    GMT_Point or None
        The site if found, else None.
    """
    for M in self.sites:
        if M.code == code:
            return M
    return None
