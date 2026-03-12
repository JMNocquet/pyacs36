"""
Shared helpers for pyacs.gts.lib.offset package.
"""

import numpy as np


def print_no_return(str):
    import sys
    sys.stdout.write(" %s" % str)
    sys.stdout.flush()


def __fmt_date(decyear):
    import pyacs.lib.astrotime
    date_cal = pyacs.lib.astrotime.decyear2datetime(decyear).strftime("%Y-%m-%d %H:%M")
    dayno, _ut = pyacs.lib.astrotime.decyear2dayno(decyear)
    return date_cal + (" doy: %03d" % dayno)


def __nargmax(data, n, threshold):
    """
    For a given 1D numpy array, returns the index of the n largest values being larger than threshold.
    return: a 2D numpy array, first column = index, second column = value, line ordered by decreasing value
    """
    D = np.copy(data)
    i = 1
    value = 9.E9
    lindex = []
    lvalue = []
    Dmin = np.min(D)
    while (i <= n) and (value >= threshold):
        arg_largest = np.argmax(D)
        value = D[arg_largest]
        if value >= threshold:
            lindex.append(arg_largest)
            lvalue.append(value)
            i = i + 1
            D[arg_largest] = Dmin
    return np.array([lindex, lvalue]).T


def get_suspected_dates(diff_data, threshold, lcomponent='NEU', verbose=False):
    """Get the list of the largest values; these are suspected offsets."""
    ldates = []
    [median_north, median_east, median_up, sn, se, su] = np.median(np.abs(diff_data), axis=0)[1:7]
    if verbose:
        print("-- median of absolute differentiated time series N %10.3lf    E %10.3lf    U %10.3lf    mm " % (median_north*1.E3, median_east*1.E3, median_up*1.E3))
    if 'N' in lcomponent:
        lindex_north = np.where(np.abs(diff_data[:, 1]) > threshold*median_north)[0].tolist()
        if verbose:
            print("-- found %d potential offsets from day-to-day difference for component North" % len(lindex_north))
        if lindex_north != []:
            lindex_north = (np.array(check_suspected_offsets(lindex_north, verbose=verbose)),)
        ldates.append(diff_data[lindex_north][:, 0].tolist())
        if verbose:
            print("-- adding %d potential offsets for component North" % len(lindex_north))
    if 'E' in lcomponent:
        lindex_east = np.where(np.abs(diff_data[:, 2]) > threshold*median_east)[0].tolist()
        if verbose:
            print("-- found %d potential offsets from day-to-day difference for component East" % len(lindex_east))
        if lindex_east != []:
            lindex_east = (np.array(check_suspected_offsets(lindex_east, verbose=verbose)),)
        ldates.append(diff_data[lindex_east][:, 0].tolist())
        if verbose:
            print("-- adding %d potential offsets for component East" % len(lindex_east))
    if 'U' in lcomponent:
        lindex_up = np.where(np.abs(diff_data[:, 3]) > threshold*median_up)[0].tolist()
        print("-- found %d potential offsets from day-to-day difference for component Up" % len(lindex_up))
        if lindex_up != []:
            lindex_up = (np.array(check_suspected_offsets(lindex_up, verbose=verbose)),)
        ldates.append(diff_data[lindex_up][:, 0].tolist())
        if verbose:
            print("-- adding %d potential offsets for component Up" % len(lindex_up))
    if not ldates:
        return []
    flat = []
    for L in ldates:
        if isinstance(L, list):
            flat.extend(L)
        else:
            flat.append(L)
    lldates = sorted(list(set(flat)))
    return lldates


def check_suspected_offsets(lindex, verbose=False):
    """
    Check that the list of suspected index does not contain two successive values.
    In this case, it is certainly an isolated outlier and the suspected offsets are removed from the list.
    """
    if verbose:
        print('-- checking successive index as indicatoin for outliers rather than offset')
    diff_lindex = np.diff(np.array(lindex))
    lbad_index = []
    for i in np.arange(len(diff_lindex)):
        if diff_lindex[i] == 1:
            lbad_index += [i, i+1]
    if verbose:
        print("-- found %d successive index" % len(lbad_index))
    lgood_index = np.delete(lindex, lbad_index)
    return lgood_index


def __ensure_list_of_list(ll):
    """Ensures ll is a list of lists. e.g. [a,b] returns [[a,b]], [[a,b]] returns [[a,b]]."""
    if not isinstance(ll, list) or len(ll) == 0:
        return [[]]
    if not isinstance(ll[0], list):
        return [ll]
    return ll
