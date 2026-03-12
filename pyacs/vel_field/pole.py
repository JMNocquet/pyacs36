"""Euler pole and pole-removal methods for Velocity_Field."""

import numpy as np
import pyacs.lib.coordinates as Coordinates
import numpy.linalg as nplinalg
import pyacs.lib.euler as euler


def calc_pole(self, plates, method='WLS', verbose=False):
    """Perform Euler pole calculation for multiple plates.

    Parameters
    ----------
    plates : dict
        Plate name -> list of site codes for each plate.
    method : str, optional
        'WLS' (weighted least-squares) or 'L1'. Default is 'WLS'.
    verbose : bool, optional
        If True, print progress. Default is False.

    Notes
    -----
    For each plate, creates: euler_stat_<plate>.dat, <plate>.gmt, <plate>.shp,
    and euler_sum.dat with Euler pole summary.
    """
    import pyacs.lib.shapefile as shapefile_mod

    H_w = {}
    for plate, lsite in sorted(plates.items()):
        print("-- processing plate: %s with %8d sites " % (plate, len(lsite)))
        if len(lsite) < 3:
            pass
        else:
            vel_plate = self.subset(lonly=lsite)
            H_w[plate] = vel_plate.pole(method=method, exp=plate, log='euler_stat_' + plate + '.dat')
            new_vel = self.substract_pole(H_w[plate][0], 'rot')
            outgmt = plate + '.gmt'
            new_vel.write(out_file=plate + '.gmt')
            shapefile_mod.psvelo_to_shapefile(outgmt, plate, verbose=verbose)

    print('-- writing euler_sum.dat')
    euler_sum_file = open('euler_sum.dat', 'w')
    euler_sum_file.write("# input file: %s\n" % self.file_name)
    euler_sum_file.write("# Euler poles\n")
    euler_sum_file.write("#-------------------------------------------------------------------------------------------------------------\n")

    for plate, lsite in sorted(plates.items()):
        [wx, wy, wz] = H_w[plate][0]
        VCV_POLE_X = H_w[plate][1]
        (llambda, phi, omega) = euler.rot2euler(wx, wy, wz)
        (max_sigma, min_sigma, azimuth, s_omega) = euler.euler_uncertainty([wx, wy, wz], VCV_POLE_X)
        euler_sum_file.write("%-10s | %8.3lf %8.3lf %6.3lf | %5.2lf %5.2lf %5.1lf %5.3lf | %11.5G %11.5G %11.5G\n" %
                             (plate, llambda, phi, omega, max_sigma, min_sigma, azimuth, s_omega, wx, wy, wz))

    euler_sum_file.write("#-------------------------------------------------------------------------------------------------------------\n")
    euler_sum_file.write("# Relative Euler poles\n")
    euler_sum_file.write("# Euler poles\n")
    euler_sum_file.write("#-------------------------------------------------------------------------------------------------------------\n")

    for plate_ref, _lsite_ref in sorted(plates.items()):
        [wx_ref, wy_ref, wz_ref] = H_w[plate_ref][0]
        VCV_POLE_X_REF = H_w[plate_ref][1]
        for plate, lsite in sorted(plates.items()):
            if plate != plate_ref:
                [wx, wy, wz] = H_w[plate][0]
                VCV_POLE_X = H_w[plate][1]
                wx = wx - wx_ref
                wy = wy - wy_ref
                wz = wz - wz_ref
                VCV_POLE_X = VCV_POLE_X + VCV_POLE_X_REF
                (llambda, phi, omega) = euler.rot2euler(wx, wy, wz)
                (max_sigma, min_sigma, azimuth, s_omega) = euler.euler_uncertainty([wx, wy, wz], VCV_POLE_X)
                euler_sum_file.write("%-5s wrt %-5s | %8.3lf %8.3lf %6.3lf | %5.2lf %5.2lf %6.1lf %5.3lf | %11.5G %11.5G %11.5G\n" %
                                     (plate, plate_ref, llambda, phi, omega, max_sigma, min_sigma, azimuth, s_omega, wx, wy, wz))

    euler_sum_file.write("#-------------------------------------------------------------------------------------------------------------\n")
    euler_sum_file.close()


def pole(self, lexclude=[], method='WLS', exp='pLate', log=None):
    """Compute Euler pole from the velocity field.

    Parameters
    ----------
    lexclude : list, optional
        Site codes to exclude from the fit.
    method : str, optional
        'WLS' (weighted least squares), 'LS', or 'Dikin' (L1). Default is 'WLS'.
    exp : str, optional
        Experiment name prefix. Default is 'pLate'.
    log : str, optional
        Log file path. If None, uses exp + '.log'.

    Returns
    -------
    tuple
        (X, VCV_POLE_X): rotation rate vector (3,) and variance-covariance (3, 3).
    """
    np.set_printoptions(precision=2, threshold=10000, linewidth=250)

    if lexclude:
        nsites = sum(1 for M in self.sites if M.code not in lexclude)
    else:
        nsites = self.nsites()

    A = np.zeros([2 * nsites, 3], float)
    B = np.zeros([2 * nsites, 1], float)
    VCV_B = np.zeros([2 * nsites, 2 * nsites], float)

    index = 0
    H_sites = {}
    for M in self.sites:
        if M.code not in lexclude:
            H_sites[M.code] = M
            R = Coordinates.mat_rot_general_to_local(np.radians(M.lon), np.radians(M.lat))
            (x, y, z) = Coordinates.geo2xyz(np.radians(M.lon), np.radians(M.lat), M.he)
            Ai = np.zeros([3, 3], float)
            Ai[0, 1] = z
            Ai[0, 2] = -y
            Ai[1, 0] = -z
            Ai[1, 2] = x
            Ai[2, 0] = y
            Ai[2, 1] = -x
            Ai = Ai / 1000.0
            RAi = np.dot(R, Ai)
            A[2 * index, 0] = RAi[0, 0]
            A[2 * index, 1] = RAi[0, 1]
            A[2 * index, 2] = RAi[0, 2]
            A[2 * index + 1, 0] = RAi[1, 0]
            A[2 * index + 1, 1] = RAi[1, 1]
            A[2 * index + 1, 2] = RAi[1, 2]
            B[2 * index, 0] = M.Ve
            B[2 * index + 1, 0] = M.Vn
            M.assign_index(index)
            VCV_B[2 * index, 2 * index] = M.SVe ** 2
            VCV_B[2 * index + 1, 2 * index + 1] = M.SVn ** 2
            cov = M.SVe * M.SVn * M.SVen
            VCV_B[2 * index + 1, 2 * index] = cov
            VCV_B[2 * index, 2 * index + 1] = cov
            index = index + 1

    P = nplinalg.inv(VCV_B)
    ATP = np.dot(A.T, P)
    N = np.dot(ATP, A)
    M_mat = np.dot(ATP, B)
    Q = nplinalg.inv(N)
    X = np.dot(Q, M_mat)
    MP = np.dot(A, X)
    RVen = (B - MP)
    VCV_POLE_X = Q * 1.E-12
    chi2 = np.asarray(np.dot(np.dot(RVen.T, P), RVen)).ravel()[0].item()
    dof = 2 * nsites - 3
    reduced_chi2 = float(np.sqrt(chi2 / float(dof)))
    (wx, wy, wz) = (X[0, 0] * 1.E-6, X[1, 0] * 1.E-6, X[2, 0] * 1.E-6)
    (llambda, phi, omega) = euler.rot2euler(wx, wy, wz)
    (max_sigma, min_sigma, azimuth, s_omega) = euler.euler_uncertainty([wx, wy, wz], VCV_POLE_X)

    flog_name = exp + '.log' if not log else log
    flog = open(flog_name, 'w')
    flog.write("\nROTATION RATE VECTOR\n")
    flog.write("Wx (rad/yr): %11.5G +- %10.5G\n" % (wx, np.sqrt(VCV_POLE_X[0, 0])))
    flog.write("Wy (rad/yr): %11.5G +- %10.5G\n" % (wy, np.sqrt(VCV_POLE_X[1, 1])))
    flog.write("Wz (rad/yr): %11.5G +- %10.5G\n" % (wz, np.sqrt(VCV_POLE_X[2, 2])))
    flog.write("ASSOCIATED VARIANCE-COVARIANCE MATRIX (rad/yr)**2\n")
    flog.write("        Wx          Wy            Wz\n")
    flog.write("---------------------------------------\n")
    flog.write("Wx |%11.5G %11.5G %11.5G\n" % (VCV_POLE_X[0, 0], VCV_POLE_X[1, 0], VCV_POLE_X[2, 0]))
    flog.write("Wy |           %11.5G %11.5G\n" % (VCV_POLE_X[1, 1], VCV_POLE_X[2, 1]))
    flog.write("Wz |                      %11.5G\n" % VCV_POLE_X[2, 2])
    flog.write("---------------------------------------\n")
    flog.write("\nEULER POLE\n")
    flog.write("longitude (dec. degree)    : %6.2lf\n" % llambda)
    flog.write("latitude  (dec. degree)    : %6.2lf\n" % phi)
    flog.write("angular velocity (deg/Myr ): %6.3lf\n" % omega)
    flog.write("ASSOCIATED ERROR ELLIPSE\n")
    flog.write("semi major axis            : %5.2lf\n" % max_sigma)
    flog.write("semi minor axis            : %5.2lf\n" % min_sigma)
    flog.write("azimuth of semi major axis : %5.1lf\n" % azimuth)
    flog.write("std(angular velocity)      : %5.3lf\n" % s_omega)
    flog.write("\nSTATISTICS\n")
    flog.write("----------\n")
    flog.write("Number of sites     = %10d\n" % nsites)
    flog.write("Chi**2              = %10.1lf\n" % chi2)
    flog.write("Reduced Chi**2      = %10.1lf\n" % reduced_chi2)
    flog.write("Deg. of. freedom    = %10d\n" % dof)
    if dof == 0:
        flog.write("A post. var. factor = No meaning (dof=0)\n")
    else:
        sigma_0 = float(np.sqrt(chi2 / dof))
        flog.write("A post. var. factor = %10.1lf\n" % sigma_0)
    flog.write("\nRESIDUALS\n")
    flog.write("----------\n")
    flog.write("site                    pred_e     pred_n         Ve         Vn       R_ve       R_vn       S_ve       S_vn      RN_ve      RN_vn\n")
    flog.write("------------------------------------------------------------------------------------------------------------------------------------\n")

    cpt_site = 0
    rms = 0.0
    wrms = 0.0
    dwrms = 0.0
    for site_code, M in sorted(H_sites.items()):
        M = H_sites[site_code]
        if M.code not in lexclude:
            i = M.get_index()
            mpe = MP[2 * i, 0]
            mpn = MP[2 * i + 1, 0]
            oe = B[2 * i, 0]
            on = B[2 * i + 1, 0]
            rve = RVen[2 * i, 0]
            rvn = RVen[2 * i + 1, 0]
            sve = M.SVe
            svn = M.SVn
            flog.write("%-19s %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf\n" %
                       (M.code, mpe, mpn, oe, on, rve, rvn, sve, svn, rve / sve, rvn / svn))
            cpt_site = cpt_site + 1
            rms = rms + rve ** 2 + rvn ** 2
            wrms = wrms + (rve / sve) ** 2 + (rvn / svn) ** 2
            dwrms = dwrms + 1.0 / sve ** 2 + 1.0 / svn ** 2
    flog.write("------------------------------------------------------------------------------------------------------------------------------------\n")
    rms = np.sqrt(rms / 2 / cpt_site)
    wrms = np.sqrt(wrms / dwrms)
    flog.write("rms = %10.2lf mm/yr    wrms = %10.2lf mm/yr\n\n" % (rms, wrms))
    flog.close()
    print("-- Rotation rate results logged in %s" % flog_name)
    return (X * 1.E-6, VCV_POLE_X)


def substract_pole(self, W=None, type_euler='rot'):
    """Subtract (or add) the prediction of an Euler pole from the velocity field.

    Parameters
    ----------
    W : array_like, optional
        Euler pole: either [lon, lat, omega] in deg/Myr or [wx, wy, wz] in rad/yr.
    type_euler : str, optional
        'rot' for cartesian rotation rate (rad/yr), 'euler' for geographic. Default is 'rot'.

    Returns
    -------
    Velocity_Field
        New velocity field with pole removed.
    """
    W = np.array(W).reshape(3, 1)
    ltype = ['euler', 'rot']
    if type_euler not in ltype:
        raise ValueError("!!! type_euler should be either euler or rot. type_euler=%s" % type_euler)

    if type_euler == 'euler':
        (lon, lat, omega) = (W[0, 0], W[1, 0], W[2, 0])
        (wx, wy, wz) = euler.euler2rot(lon, lat, omega)
        W = np.zeros([3, 1])
        W[0, 0] = wx
        W[1, 0] = wy
        W[2, 0] = wz

    lGpoint = []
    for M in self.sites:
        R = Coordinates.mat_rot_general_to_local(np.radians(M.lon), np.radians(M.lat))
        (x, y, z) = Coordinates.geo2xyz(np.radians(M.lon), np.radians(M.lat), M.he)
        Ai = np.zeros([3, 3], float)
        Ai[0, 1] = z
        Ai[0, 2] = -y
        Ai[1, 0] = -z
        Ai[1, 2] = x
        Ai[2, 0] = y
        Ai[2, 1] = -x
        RAi = np.dot(R, Ai)
        Pi = np.dot(RAi, W)
        pve = Pi[0, 0] * 1.E3
        pvn = Pi[1, 0] * 1.E3
        N = M.copy()
        N.Ve = M.Ve - pve
        N.Vn = M.Vn - pvn
        lGpoint.append(N)

    return self.__class__(lgmt_points=lGpoint)
