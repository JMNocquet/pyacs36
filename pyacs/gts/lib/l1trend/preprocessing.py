"""
Pre-processing functions for L1-trend analysis.
"""

import numpy as np
import pyacs.lib.astrotime as at

import logging
import pyacs.message.message as MESSAGE
import pyacs.message.verbose_message as VERBOSE
import pyacs.message.error as ERROR
import pyacs.message.warning as WARNING
import pyacs.message.debug_message as DEBUG


def pre_process_test(gts, component, threshold=5000., logger=None):
    """
    Pre-process. Test whether a large offset is present in the time series.
    Find the largest offset exceeding threshold found in the time series.
    Returns the index in the time series immediately after the largest offset or a void list if no offset is found

    Parameters
    ----------
    gts : Gts
        Input time series
    component : str
        Component to test ('E', 'N', or 'U')
    threshold : float
        Threshold for offset detection
    logger : logging.Logger, optional
        Logger instance for logging messages

    Returns
    -------
    list
        List of indices where offsets are found
    """
    threshold_up = threshold * 5

    # Ensure threshold is float
    threshold = float(threshold)

    # get time period
    str_sdate = at.decyear2datetime(gts.data[0,0]).isoformat()[:10]
    str_edate = at.decyear2datetime(gts.data[-1,0]).isoformat()[:10]
    str_period = f"[{str_sdate} - {str_edate}]"

    mesg = f"Running pre_process_test on Gts {gts.code} {str_period} {component} with threshold {threshold:.1f}"
    VERBOSE(mesg)
    if logger:
        logger.info(mesg)

    # case only one value
    if gts.data.shape[0] == 1:
        mesg = f"In pre_process_test on Gts {gts.code} {str_period} {component}. Only one value found. No offset returned."
        if logger:
            logger.warning(mesg)
        else:
            WARNING(mesg)
        return []

    # case only two values
    if gts.data.shape[0] == 2:
        mesg = f"In pre_process_test on Gts {gts.code} {str_period} {component}. Only two values found."
        if logger:
            logger.warning(mesg)
        else:
            WARNING(mesg)
        
        ivel = gts.ivel()
        if 'E' in component or 'N' in component:
            if np.sqrt(ivel.data[0,1]**2 + ivel.data[0,2]**2) * 1.E3 > threshold:
                str_date_offset = f"{at.decyear2datetime(gts.data[0, 0]).isoformat()[:19]}-{at.decyear2datetime(gts.data[1, 0]).isoformat()[:19]}"
                mesg = f"In pre_process_test on Gts {gts.code} {str_period} {component}. Offset detected on horizontal component {str_date_offset}"
                if logger:
                    logger.warning(mesg)
                else:
                    WARNING(mesg)
                return [1]
        if 'U' in component:
            if np.fabs(ivel.data[0,3]) * 1.E3 > threshold_up:
                str_date_offset = f"{at.decyear2datetime(gts.data[0, 0]).isoformat()[:19]}-{at.decyear2datetime(gts.data[1, 0]).isoformat()[:19]}"
                mesg = f"In pre_process_test on Gts {gts.code} {str_period} {component}. Offset detected on up component {str_date_offset}"
                if logger:
                    logger.warning(mesg)
                else:
                    WARNING(mesg)
                return [1]
        return []

    # case 3 values
    if gts.data.shape[0] == 3:
        mesg = f"In pre_process_test on Gts {gts.code} {str_period} {component}. Only three values found."
        if logger:
            logger.warning(mesg)
        else:
            WARNING(mesg)
        ivel = gts.ivel()
        if 'E' in component or 'N' in component:
            idx = np.argmax(np.sqrt(ivel.data[:,1]**2 + ivel.data[:,2]**2))
            if np.sqrt(ivel.data[idx,1]**2 + ivel.data[idx,2]**2) * 1.E3 > threshold:
                str_date_offset = f"{at.decyear2datetime(gts.data[idx, 0]).isoformat()[:19]}-{at.decyear2datetime(gts.data[idx+1, 0]).isoformat()[:19]}"
                mesg = f"In pre_process_test on Gts {gts.code} {str_period} {component}. Offset detected on horizontal component {str_date_offset}"
                if logger:
                    logger.warning(mesg)
                else:
                    WARNING(mesg)
                return [idx+1]
            else:
                mesg = f"No offset found for Gts {gts.code} {str_period} {component} horizontal component"
                VERBOSE(mesg)
                if logger:
                    logger.info(mesg)
                return []
        if 'U' in component:
            idx = np.argmax(np.fabs(ivel.data[:,3]))
            if np.fabs(ivel.data[idx,3]) * 1.E3 > threshold_up:
                str_date_offset = f"{at.decyear2datetime(gts.data[idx, 0]).isoformat()[:19]}-{at.decyear2datetime(gts.data[idx+1, 0]).isoformat()[:19]}"
                mesg = f"In pre_process_test on Gts {gts.code} {str_period} {component}. Offset detected on up component {str_date_offset}"
                if logger:
                    logger.warning(mesg)
                else:
                    WARNING(mesg)
                return [idx+1]
            else:
                mesg = f"No offset found for Gts {gts.code} {str_period} {component} horizontal component"
                VERBOSE(mesg)
                if logger:
                    logger.info(mesg)
                return []
    
    # case 4-5 values
    if gts.data.shape[0] == 4 or gts.data.shape[0] == 5:
        mesg = f"In pre_process_test on Gts {gts.code} {str_period} {component}. Only four values found."
        if logger:
            logger.warning(mesg)
        else:
            WARNING(mesg)
        its = gts.median_filter(3)
        # compute ivel
        ivel = its.ivel()
        if 'E' in component or 'N' in component:
            norm_ivel = np.sqrt(ivel.data[:, 1] ** 2 + ivel.data[:, 2] ** 2) * 1.E3
            idx = np.argmax(norm_ivel)
            if norm_ivel[idx] > threshold:
                str_date_offset = f"{at.decyear2datetime(gts.data[idx, 0]).isoformat()[:19]}-{at.decyear2datetime(gts.data[idx+1, 0]).isoformat()[:19]}"
                mesg = f"Found offset on horizontal component at middle of {str_date_offset}"
                VERBOSE(mesg)
                if logger:
                    logger.info(mesg)
                return [idx+1]
            else:
                mesg = f"No offset found for Gts {gts.code} {str_period} {component} horizontal component"
                VERBOSE(mesg)
                if logger:
                    logger.info(mesg)
                return []
        if 'U' in component:
            norm_ivel = np.fabs(ivel.data[:, 3]) * 1.E3
            idx = np.argmax(norm_ivel)
            if norm_ivel[idx] > threshold_up:
                str_date_offset = f"{at.decyear2datetime(gts.data[idx, 0]).isoformat()[:19]}-{at.decyear2datetime(gts.data[idx+1, 0]).isoformat()[:19]}"
                mesg = f"Found offset on up component at middle of {str_date_offset}"
                VERBOSE(mesg)
                if logger:
                    logger.info(mesg)
                return [idx+1]
            else:
                mesg = f"No offset found for Gts {gts.code} {str_period} {component} up component"
                VERBOSE(mesg)
                if logger:
                    logger.info(mesg)
                return []
    
    # case >= 5 values
    if gts.data.shape[0] > 5:
        its3 = gts.median_filter(3)
        its5 = gts.median_filter(5)
        # compute ivel
        ivel3 = its3.ivel()
        ivel5 = its5.ivel()
        # compute the mean of the two ivels
        if 'E' in component or 'N' in component:
            mean_ivel = np.sqrt( ((ivel3.data[:, 1] + ivel5.data[:, 1]) / 2)**2 + ((ivel3.data[:, 2] + ivel5.data[:, 2]) / 2)**2) * 1.E3
            idx = np.argmax(mean_ivel)
            if mean_ivel[idx]**2 > threshold**2:
                str_date_offset = f"{at.decyear2datetime(gts.data[idx, 0]).isoformat()[:19]}-{at.decyear2datetime(gts.data[idx+1, 0]).isoformat()[:19]}"
                mesg = f"Found offset on horizontal component at middle of {str_date_offset}"
                VERBOSE(mesg)
                if logger:
                    logger.info(mesg)
                return [idx+1]
            else:
                mesg = f"No offset found for Gts {gts.code} {str_period} {component} horizontal component"
                VERBOSE(mesg)
                if logger:
                    logger.info(mesg)
                return []
        if 'U' in component:
            mean_ivel = (ivel3.data[:, 3] + ivel5.data[:, 3]) / 2 * 1.E3
            idx = np.argmax(mean_ivel**2)
            if mean_ivel[idx]**2 > threshold_up**2:
                str_date_offset = f"{at.decyear2datetime(gts.data[idx, 0]).isoformat()[:19]}-{at.decyear2datetime(gts.data[idx+1, 0]).isoformat()[:19]}"
                mesg = f"Found offset on up component at middle of {str_date_offset}"
                VERBOSE(mesg)
                if logger:
                    logger.info(mesg)
                return [idx+1]
            else:
                mesg = f"No offset found for Gts {gts.code} {str_period} {component} up component."
                VERBOSE(mesg)
                if logger:
                    logger.info(mesg)
                return []


def pre_process_ts(self, threshold=5000.):
    """
    Pre-process a time series to avoid convergence problems in l1trend.
    It performs change point detection analysis, correct for unrealistically large velocity
    and return a cleaned time series
    """
    VERBOSE("Running pre_process_ts on Gts %s with threshold %.1lf" % (self.code, threshold))

    # detrend
    dts = self.detrend_median()
    # working time series
    nts = dts.copy()
    its = dts.median_filter(5)
    # compute ivel
    ivel = its.ivel()
    # set all reasonable ivel to zero
    # horizontal components
    lidxe = np.where(np.fabs(ivel.data[:, 1]) > threshold * 1.E-3)[0]
    lidxn = np.where(np.fabs(ivel.data[:, 2]) > threshold * 1.E-3)[0]
    lidx = np.unique(np.sort(np.append(lidxe, lidxn)))
    VERBOSE("Found %d offsets in horizontal components" % lidx.shape[0])
    for i in lidx:
        VERBOSE("%s " % (at.decyear2datetime(its.data[i, 0]).isoformat()))
    lidxx = np.zeros(ivel.data.shape[0])
    lidxx[lidx] = 1
    ivel.data[:, 1] = ivel.data[:, 1] * lidxx
    ivel.data[:, 2] = ivel.data[:, 2] * lidxx

    # vertical component
    lidx = np.where(np.fabs(ivel.data[:, 2]) > threshold * 1.E-3)[0]
    VERBOSE("Found %d offsets in up component" % lidx.shape[0])
    for i in lidx:
        VERBOSE("%s " % (at.decyear2datetime(its.data[i, 0]).isoformat()))
    lidxx = np.zeros(ivel.data.shape[0])
    lidxx[lidx] = 1
    ivel.data[:, 3] = ivel.data[:, 3] * lidxx

    # compute its
    for i in [1, 2, 3]:
        its.data[:, i] = np.append(0, np.cumsum(ivel.data[:, i] * np.diff(its.data[:, 0])))
        nts.data[:, i] = dts.data[:, i] - its.data[:, i]
    # remove obvious outliers
    ots = nts.find_outliers_simple()
    its.outliers = ots.outliers
    # return
    return ots.remove_outliers(), its.remove_outliers()
