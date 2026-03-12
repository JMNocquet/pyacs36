"""Strain rate methods for Velocity_Field."""

import numpy as np
import pyacs.lib.glinalg as glinalg


def vgrad(X, Y, VE, VN, SVE=None, SVN=None, CVEN=None, CODE=None, method='WLS', verbose=False):
    """Linear estimate of horizontal velocity gradient from site velocities.

    Parameters
    ----------
    X, Y : numpy.ndarray
        1D arrays of longitudes and latitudes (decimal degrees).
    VE, VN : numpy.ndarray
        East and north velocities (mm/yr).
    SVE, SVN : numpy.ndarray, optional
        Standard deviations of VE, VN (mm/yr).
    CVEN : numpy.ndarray, optional
        Correlation or covariance term for VE/Vn.
    CODE : numpy.ndarray, optional
        Site codes (strings).
    method : str, optional
        'WLS' (weighted least-squares) or 'L1'. Default is 'WLS'.
    verbose : bool, optional
        If True, print progress. Default is False.

    Returns
    -------
    tuple
        l0, p0, SOL, COV, RESIDUALS, chi2: barycenter (lon, lat), solution vector,
        covariance, residuals, chi-square.
    """
    import pyacs.lib.coordinates

    def __vgrad_obs_eq__(l0, p0, l, p, ve, vn, sve, svn, cven):
        """Observation equation for horizontal velocity gradient.
        Observations in mm/yr; result in nstrain/yr.
        """
        Rt = 6371.E3
        Rt3nstrain = Rt / 1.E6
        Ai = np.zeros((2, 6))
        Bi = np.zeros(2)
        Corri = np.zeros((2, 2))
        Ai[0, 0] = 1.
        Ai[1, 1] = 1.
        Ai[0, 2] = np.radians(l - l0) * Rt3nstrain * np.cos(np.radians((p + p0) / 2.))
        Ai[0, 3] = np.radians(p - p0) * Rt3nstrain
        Ai[1, 4] = np.radians(l - l0) * Rt3nstrain * np.cos(np.radians((p + p0) / 2.))
        Ai[1, 5] = np.radians(p - p0) * Rt3nstrain
        Bi[0] = ve
        Bi[1] = vn
        Corri[0, 0] = 1.
        Corri[1, 1] = 1.
        Corri[1, 0] = cven
        Corri[0, 1] = cven
        Cvi = glinalg.corr_to_cov(Corri, np.array([sve, svn]))
        return (Ai, Bi, Cvi)

    def __barycenter__(X_, Y_):
        """Return the barycenter of points in geographical coordinates."""
        XYZ = np.array(list(map(pyacs.lib.coordinates.geo2xyz, np.radians(X_), np.radians(Y_), np.zeros(X_.size))))
        [mx, my, mz] = np.mean(XYZ, axis=0)
        [mlong, mlat, _] = pyacs.lib.coordinates.xyz2geo(mx, my, mz, unit='dec_deg')
        return (mlong, mlat)

    (l0, p0) = __barycenter__(X, Y)
    if verbose:
        print("-- Barycenter at (%10.5lf , %10.5lf )" % (l0, p0))

    A = np.zeros((2 * X.size, 6))
    B = np.zeros((2 * X.size))
    CV = np.zeros((2 * X.size, 2 * X.size))

    for i in np.arange(X.size):
        (Ai, Bi, Cvi) = __vgrad_obs_eq__(l0, p0, X[i], Y[i], VE[i], VN[i], SVE[i], SVN[i], CVEN[i])
        A[2 * i:2 * i + 2, :] = Ai
        B[2 * i:2 * i + 2] = Bi
        CV[2 * i:2 * i + 2, 2 * i:2 * i + 2] = Cvi

    SOL, COV, RESIDUALS, chi2 = glinalg.lscov_full(A, B, CV)
    return (l0, p0, SOL, COV, RESIDUALS, chi2)


def strain(self, lcode, save=None, method='WLS', verbose=False):
    """Calculate strain rate from a list of sites.

    Parameters
    ----------
    lcode : list
        List of site codes to use.
    save : str, optional
        File path to save the result. If None, result is only printed.
    method : str, optional
        'WLS' (weighted least-squares) or 'L1'. Default is 'WLS'.
    verbose : bool, optional
        If True, print progress. Default is False.

    Returns
    -------
    Velocity_Field
        self.
    """
    X = np.array([])
    Y = np.array([])
    VE = np.array([])
    VN = np.array([])
    SVE = np.array([])
    SVN = np.array([])
    CVEN = np.array([])
    CODE = np.array([])

    for code in lcode:
        M = self.site(code)
        X = np.append(X, M.lon)
        Y = np.append(Y, M.lat)
        VE = np.append(VE, M.Ve)
        VN = np.append(VN, M.Vn)
        SVE = np.append(SVE, M.SVe)
        SVN = np.append(SVN, M.SVn)
        CVEN = np.append(CVEN, M.SVen)
        CODE = np.append(CODE, code)

    l0, p0, SOL, COV, RESIDUALS, chi2 = vgrad(X, Y, VE, VN, SVE, SVN, CVEN, CODE, method='WLS', verbose=verbose)
    CORR, SIGMA = glinalg.cov_to_corr(COV)

    [ve0, vn0, dvx_dx, dvx_dy, dvy_dx, dvy_dy] = SOL
    [sve0, svn0, sdvx_dx, sdvx_dy, sdvy_dx, sdvy_dy] = SIGMA
    sven0 = CORR[0, 1]

    TRANS_GRADV = np.zeros((6, 6))
    TRANS_GRADV[0, 2] = 1.
    TRANS_GRADV[1, 3] = .5
    TRANS_GRADV[1, 4] = .5
    TRANS_GRADV[2, 5] = 1.
    TRANS_GRADV[3, 3] = .5
    TRANS_GRADV[3, 4] = -.5
    TRANS_GRADV[4, 2] = 1.
    TRANS_GRADV[4, 5] = -1.
    TRANS_GRADV[5, 3] = 1.
    TRANS_GRADV[5, 4] = 1.

    TRANS = np.dot(TRANS_GRADV, SOL)
    [edot11, edot12, edot22, omega, gamma1, gamma2] = TRANS

    VAR_TRANS = np.dot(TRANS_GRADV, np.dot(COV, TRANS_GRADV.T))
    COV_STRAIN_RATE = VAR_TRANS[:3, :3]

    s_omega = np.sqrt(VAR_TRANS[3, 3])
    s_gamma1 = np.sqrt(VAR_TRANS[4, 4])
    s_gamma2 = np.sqrt(VAR_TRANS[5, 5])

    eps1 = .5 * (SOL[2] + SOL[5] + np.sqrt((SOL[2] - SOL[5]) ** 2 + (SOL[3] + SOL[4]) ** 2))
    eps2 = .5 * (SOL[2] + SOL[5] - np.sqrt((SOL[2] - SOL[5]) ** 2 + (SOL[3] + SOL[4]) ** 2))
    azimut = .5 * np.degrees(np.arctan((SOL[3] + SOL[4]) / (SOL[5] - SOL[2]))) + 90.0

    gamma = np.sqrt(gamma1 ** 2 + gamma2 ** 2)

    term = .5 / np.sqrt((SOL[2] - SOL[5]) ** 2 + (SOL[3] + SOL[4]) ** 2)
    dgdl11 = term * 2. * (SOL[2] - SOL[5])
    dgdl12 = term * 2. * (SOL[3] + SOL[4])
    dgdl21 = term * 2. * (SOL[3] + SOL[4])
    dgdl22 = -term * 2. * (SOL[2] - SOL[5])

    bigg = np.zeros((4, 6))
    bigg[0, 2] = +1. * (.5 + .5 * dgdl11)
    bigg[0, 3] = +1. * (.5 * dgdl12)
    bigg[0, 4] = +1. * (.5 * dgdl21)
    bigg[0, 5] = +1. * (.5 + .5 * dgdl22)
    bigg[1, 2] = +1. * (.5 - .5 * dgdl11)
    bigg[1, 3] = -1. * (.5 * dgdl12)
    bigg[1, 4] = -1. * (.5 * dgdl21)
    bigg[1, 5] = +1. * (.5 - .5 * dgdl22)
    term1 = (SOL[3] + SOL[4]) / (SOL[5] - SOL[2])
    term2 = 1. + term1 ** 2
    term3 = SOL[5] - SOL[2]
    bigg[2, 2] = 0.5 * term1 / term3 / term2
    bigg[2, 3] = 0.5 / term3 / term2
    bigg[2, 4] = 0.5 / term3 / term2
    bigg[2, 5] = -0.5 * term1 / term3 / term2
    bigg[3, 2] = dgdl11
    bigg[3, 3] = dgdl12
    bigg[3, 4] = dgdl21
    bigg[3, 5] = dgdl22

    UNC = np.dot(bigg, np.dot(COV, bigg.T))
    [seps1, seps2, sazimut, sgamma] = np.sqrt(np.diag(UNC))

    dof = 2 * X.size - 6
    if dof != 0:
        var_post = np.sqrt(chi2 / dof)
    else:
        var_post = 0.0

    wrms = np.sqrt(chi2 / (np.sum(1. / SVE) + np.sum(1. / SVN)))

    SRV = np.array([edot11, edot12, edot22]).reshape(-1, 1)
    chi2_strain_rate = np.dot(SRV.T, np.dot(glinalg.cov_to_invcov(COV_STRAIN_RATE), SRV))
    chi2_rotation_rate = omega ** 2 / s_omega ** 2

    RS = ("INFORMATION\n")
    RS += ("----------\n")
    RS += ("file = %s\n" % self.file_name)
    RS += ("Number of sites used in strain rate calculation: %d \n" % (X.size))
    RS += ("%s \n" % (' ').join(lcode))
    RS += ("Velocity at barycenter %10.5lf %10.5lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf\n" % (l0, p0, ve0, vn0, sve0, svn0, sven0))
    RS += ("\n")
    RS += ("VELOCITY GRADIENT TENSOR VALUES\n")
    RS += ("-------------------------------\n")
    RS += ("  dVx/dx     dVx/dy     dVy/dx      dVy/dy\n")
    RS += ("----------  ---------  ----------  ---------\n")
    RS += ("%10.3lf %10.3lf %10.3lf %10.3lf\n" % (dvx_dx, dvx_dy, dvy_dx, dvy_dy))
    RS += ("VELOCITY GRADIENT TENSOR UNCERTAINTIES\n")
    RS += ("--------------------------------------\n")
    RS += (" sdVx/dx     sdVx/dy    sdVy/dx     sdVy/dy\n")
    RS += ("----------   --------   ---------    --------\n")
    RS += ("%10.3lf %10.3lf %10.3lf %10.3lf\n" % (sdvx_dx, sdvx_dy, sdvy_dx, sdvy_dy))
    RS += ("Unit: nstrain / year \n")
    RS += ("FULL COVARIANCE MATRIX (not rescaled)\n")
    RS += ("-------------------------------------\n")
    RS += ("ve     %10.5lf\n" % (COV[0, 0]))
    RS += ("vn     %10.5lf %10.5lf\n" % tuple(COV[1, :2]))
    RS += ("dvx/dx %10.5lf %10.5lf %10.5lf\n" % tuple(COV[2, :3]))
    RS += ("dvx/dy %10.5lf %10.5lf %10.5lf %10.5lf\n" % tuple(COV[3, :4]))
    RS += ("dvy/dx %10.5lf %10.5lf %10.5lf %10.5lf %10.5lf \n" % tuple(COV[4, :5]))
    RS += ("dvy/dy %10.5lf %10.5lf %10.5lf %10.5lf %10.5lf %10.5lf\n" % tuple(COV[5, :6]))
    RS += ("\n")
    RS += ("STRAIN RATE TENSOR\n")
    RS += ("------------------\n")
    RS += ("     Exx         Exy       Eyy\n")
    RS += ("  --------   --------  ---------\n")
    RS += ("%10.3lf %10.3lf %10.3lf \n" % tuple(TRANS[:3]))
    RS += ("STRAIN RATE COVARIANCE\n")
    RS += ("----------------------\n")
    RS += ("          Exx          Exy           Eyy\n")
    RS += ("Exx %10.5lf \n" % (COV_STRAIN_RATE[0, 0]))
    RS += ("Exy %10.5lf %10.5lf \n" % tuple(COV_STRAIN_RATE[1, :2]))
    RS += ("Eyy %10.5lf %10.5lf %10.5lf \n" % tuple(COV_STRAIN_RATE[2, :3]))
    RS += ("\n")
    RS += ("PRINCIPAL AXIS \n")
    RS += ("---------------\n")
    RS += ("\n")
    RS += ("Eps1 = %10.5lf +- %10.5lf \n" % (eps1, seps1))
    RS += ("Eps2 = %10.5lf +- %10.5lf \n" % (eps2, seps2))
    RS += ("Azimut : %10.5lf +- %10.5lf \n" % (azimut, np.degrees(sazimut)))
    RS += ("\n")
    RS += ("Unit: nstrain / year and decimal degree\n")
    RS += ("Eps1: most extensional eigenvalue of strain tensor\n")
    RS += ("Eps2: most compressional eigenvalue of strain tensor\n")
    RS += ("Extension is taken positive\n")
    RS += ("azimut is the one of eps2 in degrees CW from North.\n")
    RS += ("\n")
    RS += ("ROTATION\n")
    RS += ("--------\n")
    RS += ("rotation in deg. / Myr :  %10.3lf +- %10.3lf \n" % (omega * 1.E-9 * 180. / np.pi * 1.E6, s_omega * 1.E-9 * 180. / np.pi * 1.E6))
    RS += ("rotation in rad /  yr  :  %.4E +- %.4E \n" % (omega * 1.E-9, s_omega * 1.E-9))
    RS += ("\n")
    RS += ("SHEAR RATES\n")
    RS += ("-----------\n")
    RS += ("gamma1 in nstrain / yr : %10.5lf +- %10.5lf \n" % (gamma1, s_gamma1))
    RS += ("gamma2 in nstrain / yr : %10.5lf +- %10.5lf \n" % (gamma2, s_gamma2))
    RS += ("gamma  in nstrain / yr : %10.5lf +- %10.5lf \n" % (gamma, sgamma))
    RS += ("\n")
    RS += ("STATISTICS\n")
    RS += ("----------\n")
    RS += ("chi2                       : %10.5lf \n" % (chi2))
    RS += ("degree of freedom          : %10d \n" % (dof))
    RS += ("posterior variance factor  : %10.5lf \n" % (var_post))
    RS += ("wrms (mm/yr)               : %10.5lf \n" % (wrms))
    RS += ("\n")
    RS += ("TEST OF SIGNIFICANCY\n")
    RS += ("--------------------\n")
    RS += ("\n")
    RS += ("tested              chi2      threshold values\n")
    RS += ("parameter           values       95%       99%\n")
    RS += ("-----------------------------------------------\n")
    RS += ("strain rate tensor %10.2lf   7.82     11.35\n" % float(np.asarray(chi2_strain_rate).flat[0]))
    RS += ("rotation rate      %10.2lf   3.84      6.64\n" % float(np.asarray(chi2_rotation_rate).flat[0]))
    RS += ("\n")
    RS += ("RESIDUALS\n")
    RS += ("---------\n")
    RS += ("code      long.       lat.       Ve        Vn    Res_Ve    Res_Vn   sigma_Ve   sigma_Vn  RN_Ve     RN_Vn\n")
    RS += ("---------------------------------------------------------------------------------------------------------\n")
    for i in np.arange(X.size):
        RS += ("%s %10.5lf %10.5lf %8.2lf  %8.2lf  %8.2lf  %8.2lf  %8.2lf  %8.2lf  %8.2lf  %8.2lf  \n" %
               (CODE[i], X[i], Y[i], VE[i], VN[i], RESIDUALS[2 * i], RESIDUALS[2 * i + 1], SVE[i], SVN[i], RESIDUALS[2 * i] / SVE[i], RESIDUALS[2 * i + 1] / SVN[i]))

    if verbose or (save is None):
        print(RS)

    if save is not None:
        f = open(save, 'w')
        f.write(RS)
        f.close()

    return self
