"""
Simple empirical procedure to find offsets (multi-step: large, intermediate, small).
"""

from ._utils import __fmt_date, get_suspected_dates


def _print_dates_4digits(L):
    if isinstance(L, float):
        L = [L]
    return ["%9.4lf" % date for date in sorted(L)]


def find_offsets(self, threshold=3, n_max_offsets=9, conf_level=95, lcomponent='NE', verbose=True, in_place=False):
    """Simple empirical procedure to find offsets.

    Parameters
    ----------
    threshold : float, optional
        Threshold for preliminary offset detection. Default is 3.
    n_max_offsets : int, optional
        Maximum number of offsets to detect. Default is 9.
    conf_level : float, optional
        Confidence level (percent) to accept offset. Default is 95.
    lcomponent : str, optional
        Components for detection ('N','E','U'). Default is 'NE'.
    verbose : bool, optional
        Verbose mode. Default is True.
    in_place : bool, optional
        If True, modify self. Default is False.

    Returns
    -------
    Gts
        New Gts (or self if in_place) with offsets_dates and outliers set.
    """
    tts = self.copy()
    loutliers_dates = []
    gross_threshold = 10

    if verbose:
        print("********************************************************************")
        print("-- %s STEP #1: trying to identify large offsets" % self.code)
        print("********************************************************************")

    tmp_gts_gross = tts.suspect_offsets(threshold=gross_threshold, verbose=True, lcomponent=lcomponent, n_max_offsets=n_max_offsets)
    if verbose:
        print("-- %1d suspected large offsets found" % (len(tmp_gts_gross.offsets_dates)))
        for odate in tmp_gts_gross.offsets_dates:
            print("-- %s %12.8lf" % (__fmt_date(odate), odate))

    significant_offsets = []
    for offset_date in tmp_gts_gross.offsets_dates:
        if tts.test_offset_significance(offset_date, conf_level=conf_level, lcomponent=lcomponent, verbose=verbose, mode='local') \
           and tts.test_offset_significance(offset_date, conf_level=conf_level, lcomponent=lcomponent, verbose=verbose, mode='detrend'):
            significant_offsets.append(offset_date)
        else:
            if verbose:
                print("=> Since it is not an offset, date %10.4lf is potentially an outlier" % offset_date)
            tmp_gts_gross.find_outlier_around_date(offset_date, conf_level=conf_level, n=3, lcomponent=lcomponent, verbose=verbose)

    if verbose:
        print("=> Large offset search : %1d offsets confirmed, %1d were actually outliers" % (len(significant_offsets), len(tmp_gts_gross.outliers)))

    tts.offsets_dates += significant_offsets
    if tmp_gts_gross.outliers != []:
        lindex_outliers = tmp_gts_gross.outliers
        new_outliers_dates = tmp_gts_gross.data[lindex_outliers, 0].tolist()
        loutliers_dates = new_outliers_dates

    loutliers_gross = tmp_gts_gross.find_outliers_simple(threshold=gross_threshold).outliers
    loutliers_dates += self.data[loutliers_gross, 0].tolist()
    if verbose:
        print("=> Outliers search : %03d outliers found" % (len(loutliers_dates)))

    tmp_gts_gross.plot()

    if verbose:
        print("********************************************************************")
        print("-- %s STEP #2: trying to identify intermediate size offsets" % self.code)
        print("********************************************************************")

    intermediate_threshold = 0.5 * (10. + threshold)
    previous_offsets_values = tmp_gts_gross.remove_outliers().detrend().offsets_values
    tmp_gts_intermediate = tmp_gts_gross.remove_outliers().apply_offsets(previous_offsets_values).suspect_offsets(threshold=intermediate_threshold, verbose=True, lcomponent='NE', n_max_offsets=10)
    str_offsets_dates = " ".join(_print_dates_4digits(tmp_gts_intermediate.offsets_dates))
    if verbose:
        print("=> %1d potential intermediate offsets found at %s" % (len(tmp_gts_intermediate.offsets_dates), str_offsets_dates))

    significant_offsets = []
    lpotential_offsets_dates = tmp_gts_intermediate.offsets_dates
    tmp_gts_intermediate.offsets_dates = []

    for offset_date in lpotential_offsets_dates:
        if tmp_gts_intermediate.test_offset_significance(offset_date, conf_level=conf_level, lcomponent=lcomponent, verbose=verbose, mode='local') \
           and tmp_gts_intermediate.test_offset_significance(offset_date, conf_level=conf_level, lcomponent=lcomponent, verbose=verbose, mode='detrend'):
            significant_offsets.append(offset_date)
        else:
            if verbose:
                print("=> Since it is not an offset, date %10.4lf is potentially an outlier" % offset_date)
            tmp_gts_intermediate.find_outlier_around_date(offset_date, conf_level=conf_level, n=3, lcomponent=lcomponent, verbose=verbose)

    if verbose:
        print("=> Intermediate size offsets search : %1d offsets confirmed, %1d were actually outliers" % (len(significant_offsets), len(tmp_gts_intermediate.outliers)))

    tts.offsets_dates += significant_offsets
    if tmp_gts_intermediate.outliers != []:
        lindex_outliers = tmp_gts_intermediate.outliers
        new_outliers_dates = tmp_gts_intermediate.data[lindex_outliers, 0].tolist()
        loutliers_dates = loutliers_dates + new_outliers_dates

    loutliers_intermediate = tmp_gts_intermediate.find_outliers_simple(threshold=intermediate_threshold).outliers
    loutliers_dates += self.data[loutliers_intermediate, 0].tolist()
    if verbose:
        print("=> Outliers search : %03d outliers found" % (len(loutliers_dates)))

    if verbose:
        print("********************************************************************")
        print("-- %s STEP #3: trying to identify small offsets" % self.code)
        print("********************************************************************")

    previous_offsets_values = tmp_gts_gross.remove_outliers().detrend().offsets_values
    tmp_gts_final = tmp_gts_intermediate.remove_outliers().apply_offsets(previous_offsets_values).suspect_offsets(threshold=threshold, verbose=True, lcomponent='NE', n_max_offsets=10)
    str_offsets_dates = " ".join(_print_dates_4digits(tmp_gts_final.offsets_dates))
    if verbose:
        print("=> %1d potential subtle offsets found at %s" % (len(tmp_gts_final.offsets_dates), str_offsets_dates))

    significant_offsets = []
    lpotential_offsets_dates = tmp_gts_final.offsets_dates
    tmp_gts_final.offsets_dates = []

    for offset_date in lpotential_offsets_dates:
        if tmp_gts_intermediate.test_offset_significance(offset_date, conf_level=conf_level, lcomponent=lcomponent, verbose=verbose, mode='local') \
           and tmp_gts_intermediate.test_offset_significance(offset_date, conf_level=conf_level, lcomponent=lcomponent, verbose=verbose, mode='detrend_seasonal'):
            significant_offsets.append(offset_date)
        else:
            tmp_gts_final.find_outlier_around_date(offset_date, conf_level=conf_level, n=3, lcomponent=lcomponent, verbose=verbose)

    if verbose:
        print("=> Subtle offset search: %1d offsets confirmed, %1d were actually outliers" % (len(significant_offsets), len(tmp_gts_final.outliers)))

    tts.offsets_dates += significant_offsets
    offsets_dates = sorted(list(set(tts.offsets_dates)))
    tts.offsets_dates = offsets_dates

    if tmp_gts_final.outliers != []:
        lindex_outliers = tmp_gts_final.outliers
        new_outliers_dates = tmp_gts_final.data[lindex_outliers, 0].tolist()
        loutliers_dates = loutliers_dates + new_outliers_dates

    loutliers_final = tmp_gts_final.find_outliers_simple(threshold=threshold).outliers
    loutliers_dates += self.data[loutliers_final, 0].tolist()
    if verbose:
        print("=> Outliers search : %03d outliers found" % (len(loutliers_dates)))

    if loutliers_dates != []:
        from pyacs.gts.Gts import get_index_from_dates
        returned_index = get_index_from_dates(loutliers_dates, tts.data, tol=0.25)
        loutliers = sorted(list(set(returned_index)))
    else:
        loutliers = []

    if verbose:
        print("**********************************************************************************")
        print("=> Final results of offset search: %02d offsets found; additionally %02d outliers were flagged" % (len(offsets_dates), len(loutliers_dates)))
        print("**********************************************************************************")

    new_Gts = self.copy()
    new_Gts.offsets_dates = offsets_dates
    new_Gts.outliers = loutliers

    if in_place:
        self.offsets_dates = new_Gts.offsets_dates
        self.outliers = new_Gts.outliers
        return self
    else:
        return new_Gts
