"""
Test whether an offset at a given date is statistically significant.
"""

import numpy as np

from ._utils import __fmt_date, print_no_return


def test_offset_significance(self, date, conf_level=95, lcomponent='NE', verbose=True, debug=False, mode='local'):
    """Test whether an offset is statistically significant.

    Parameters
    ----------
    date : float
        Offset date in decimal year.
    conf_level : float, optional
        Confidence level in percent. Default is 95.
    lcomponent : str, optional
        Components to test ('N','E','U'). Default is 'NE'.
    verbose : bool, optional
        Verbose mode. Default is True.
    debug : bool, optional
        Debug output. Default is False.
    mode : str, optional
        'local', 'detrend', or 'detrend_seasonal'. Default is 'local'.

    Returns
    -------
    bool
        True if significant, else False.
    """
    def f_ratio(chi_square_1, p1, chi_square_2, p2, n):
        F = ((chi_square_1 - chi_square_2) / (p2 - p1)) / (chi_square_2 / (n - p2))
        from scipy.stats import f
        return f.cdf(F, p2 - p1, n - p2)

    if verbose:
        print("-- Testing offset for site %s date %s using mode %s on components %s with confidence level %6.1lf %%" % (self.code, __fmt_date(date), mode, lcomponent, conf_level))

    OK = False
    component = {}
    component[1] = 'North'
    component[2] = 'East'
    component[3] = 'Up'
    li = []
    if 'N' in lcomponent:
        li.append(1)
    if 'E' in lcomponent:
        li.append(2)
    if 'U' in lcomponent:
        li.append(3)
    H = {}
    score = {}

    if mode == 'smooth':
        l_n_days = [5, 7, 9, 11, 13, 15, 17]
        lwindow_len = [3, 3, 5, 5, 5, 7, 7]
        coeff = [0.5, 1, 3, 4, 10, 8, 6]
        for i in sorted(li):
            print_no_return("-- n_samples/probability %5s:" % component[i])
            k = -1
            for n in l_n_days:
                k = k + 1
                before = np.array([])
                after = np.array([])
                around = np.array([])
                try:
                    before = self.extract_ndates_before_date(date, n).smooth(window_len=lwindow_len[k]).data[:, i] - self.extract_ndates_before_date(date, n).data[:, i]
                    after = self.extract_ndates_after_date(date, n).smooth(window_len=lwindow_len[k]).data[:, i] - self.extract_ndates_after_date(date, n).data[:, i]
                    around = self.extract_ndates_around_date(date, n).smooth(window_len=lwindow_len[k]).data[:, i] - self.extract_ndates_around_date(date, n).data[:, i]
                except Exception:
                    if debug:
                        print("!!! Could not test offset significancy for date %lf using n=%d surrounding data" % (date, n))
                    break
                if debug:
                    print("-- Testing significancy for date %lf using n=%d surrounding data" % (date, n))
                if before.shape[0] * after.shape[0] * around.shape[0] != 0.0:
                    before_std = np.std(before) * 1.E3
                    after_std = np.std(after) * 1.E3
                    around_std = np.std(around) * 1.E3
                    chi_square_2 = before.shape[0] * before_std**2 + after.shape[0] * after_std**2
                    chi_square_1 = around.shape[0] * around_std**2
                    H[i, k] = f_ratio(chi_square_1, 1, chi_square_2, 2, around.shape[0]) * 100.0
                    print_no_return(" %02d %5.2lf%% " % (n, H[i, k]))
            print("")
        for i in sorted(li):
            summ = 0.0
            score[i] = 0.0
            if H != {}:
                k = -1
                for n in l_n_days:
                    k = k + 1
                    try:
                        score[i] += H[i, k] * coeff[k]
                        summ = summ + coeff[k]
                    except Exception:
                        if debug:
                            print("-- No data for n=%d " % n)
                        pass
                score[i] = score[i] / summ
            if verbose:
                print("-- Probability of %s offset for site %s date %s component %5s using model: %6.1lf %%" % (mode.upper(), self.code, __fmt_date(date), component[i], score[i]))

    if mode == 'local':
        l_n_days = range(2, 10)
        coeff = [0.5, 1, 1, 4, 10, 8, 6]
        for i in sorted(li):
            print_no_return("-- n_samples/probability %5s:" % component[i])
            for n in l_n_days:
                before = np.array([])
                after = np.array([])
                around = np.array([])
                try:
                    before = self.extract_ndates_before_date(date, n).data[:, i]
                    after = self.extract_ndates_after_date(date, n).data[:, i]
                    around = self.extract_ndates_around_date(date, n).data[:, i]
                except Exception:
                    if debug:
                        print("!!! Could not test offset significancy for date %lf using n=%d surrounding data" % (date, n))
                    break
                if debug:
                    print("-- Testing significancy for date %lf using n=%d surrounding data" % (date, n))
                if before.shape[0] * after.shape[0] * around.shape[0] != 0.0:
                    before_std = np.std(before) * 1.E3
                    after_std = np.std(after) * 1.E3
                    around_std = np.std(around) * 1.E3
                    chi_square_2 = before.shape[0] * before_std**2 + after.shape[0] * after_std**2
                    chi_square_1 = around.shape[0] * around_std**2
                    H[i, n] = f_ratio(chi_square_1, 1, chi_square_2, 2, around.shape[0]) * 100.0
                    print_no_return(" %02d %5.2lf%% " % (n, H[i, n]))
            print("")
        for i in sorted(li):
            summ = 0.0
            score[i] = 0.0
            for n in range(2, 7):
                try:
                    score[i] += H[i, n] * coeff[n]
                    summ = summ + coeff[n]
                except Exception:
                    if debug:
                        print("-- No data for n=%d " % n)
                    pass
            score[i] = score[i] / summ
            if verbose:
                print("-- Probability of %s offset for site %s date %s component %5s using model: %6.1lf %%" % (mode.upper(), self.code, __fmt_date(date), component[i], score[i]))

    if mode == 'detrend':
        for i in sorted(li):
            tmp_ts = self.copy()
            tmp_ts.offsets_dates = []
            chi_square_1 = np.std(tmp_ts.detrend().data[:, i])**2 * tmp_ts.data.shape[0]
            chi_square_2 = np.std(tmp_ts.add_offsets_dates([date]).detrend().data[:, i])**2 * tmp_ts.data.shape[0]
            score[i] = f_ratio(chi_square_1, 2, chi_square_2, 3, tmp_ts.data.shape[0]) * 100.0
            if verbose:
                print("-- Probability of %s offset for site %s date %s component %5s using model: %6.1lf %%" % (mode.upper(), self.code, __fmt_date(date), component[i], score[i]))

    if mode == 'detrend_seasonal':
        for i in sorted(li):
            chi_square_1 = np.sum(self.detrend().data[:, i]**2)
            chi_square_2 = np.sum(self.add_offsets_dates([date]).detrend_seasonal().data[:, i]**2)
            n = self.data.shape[0]
            score[i] = f_ratio(chi_square_1, 1, chi_square_2, 3, 2 * n) * 100.0
            if verbose:
                print("-- Probability of %s offset for site %s date %s component %5s using model: %6.1lf %%" % (mode.upper(), self.code, __fmt_date(date), component[i], score[i]))

    score_max = np.max(np.array(list(score.values())))
    if score_max > conf_level:
        OK = True
    if verbose:
        if OK:
            print("-- Offset IS significant at the %5.2lf%% confidence level (threshold=%5.2lf%%) " % (score_max, conf_level))
        else:
            print("-- Offset NOT significant %5.2lf%% confidence level (threshold=%5.2lf%%) " % (score_max, conf_level))
    return OK
