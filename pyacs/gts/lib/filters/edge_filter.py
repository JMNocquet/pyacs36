def edge_filter(self, alpha='auto',ndays_mf=1):
    """L1 edge filter.

    Uses https://pypi.org/project/trendfilter/

    Parameters
    ----------
    alpha : str or float, optional
        Filter parameter; 'auto' for BIC-optimal. Default is 'auto'.
    ndays_mf : int, optional
        Median filter window (days) applied before L1 edge filter. Default is 1 (no median filter).

    Returns
    -------
    Gts
        L1 edge filter as a new Gts instance.
    """

    try:
        from trendfilter import trend_filter
    except ModuleNotFoundError:
        from pathlib import Path
        theme_path = Path(__file__).resolve().parents[3] / 'ts' / 'bokeh_theme.yaml'
        if not theme_path.exists():
            theme_path.write_text(
                "attrs:\n"
                "    Axis:\n"
                "        major_label_text_font_size: '12pt'\n"
                "        axis_label_text_font_size: '14pt'\n"
                "    Line:\n"
                "        line_width: 3\n"
                "    Title:\n"
                "        text_font_size: '12pt'\n"
                "    Legend:\n"
                "        background_fill_alpha: 0.8\n",
                encoding='utf-8'
            )
        from trendfilter import trend_filter
    import numpy as np

    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG
    import pyacs.debug

    from scipy.stats import chi2

    # median filter if user required
    x = self.data[:, 0]
    if ndays_mf > 1:
        gts = self.median_filter(ndays_mf)
    else:
        gts = self.copy()

    # subroutine for optimal filtering parameter search
    def f(x, y):

        from scipy.stats import chi2

        def get_stats_l1_model(y, fy, alpha):
            n = y.shape[0]
            cchi2 = np.sum((fy - y) ** 2)
            dy = np.diff(fy) / (np.diff(x))
            cp = np.where(np.fabs(np.diff(dy)) > 0.001)[0].shape[0]
            AIC = n * np.log(cchi2 / n) + 2 * cp
            BIC = n * np.log(cchi2 / n) + np.log(n) * cp
            cdof = x.shape[0] - cp
            cchi2_probability = chi2.cdf(cchi2, cdof) * 100.
            VERBOSE("alpha: %f chi2 : %f fit %.2lf n_param %d chi2_proba: %.2lf AIC: %.2lf BIC: %.2lf " % (
            alpha, cchi2, np.sqrt(cchi2 / x.shape[0]), cp, cchi2_probability, AIC, BIC))
            if pyacs.debug():
                import matplotlib.pyplot as plt
                plt.plot(x, y, 'bo', markersize=2.)
                plt.plot(x, fy, label=("%.3lf - %.1lf - %.1lf" % (alpha, cchi2, cchi2_probability)))
            return alpha, cchi2, np.sqrt(cchi2 / x.shape[0]), cp, cchi2_probability, AIC, BIC

        # search for order of magnitude
        alpha_ref = 1.
        alpha = 1.
        fy = trend_filter(x, y, l_norm=1, alpha_1=alpha_ref, constrain_zero=False)['y_fit']
        BIC_REF = get_stats_l1_model(y, fy, alpha_ref)[-1]
        GO = True
        delta = 10
        while GO:
            alpha = alpha * delta
            fy = trend_filter(x, y, l_norm=1, alpha_1=alpha, constrain_zero=False)['y_fit']
            BIC = get_stats_l1_model(y, fy, alpha)[-1]
            if BIC > BIC_REF:
                GO = False
            else:
                BIC_REF = BIC
                alpha_ref = alpha
        VERBOSE("alpha_ref: %f BIC_REF: %f" % (alpha_ref, BIC_REF))
        GO = True
        delta = 0.1
        alpha = alpha_ref
        while GO:
            alpha = alpha * delta
            fy = trend_filter(x, y, l_norm=1, alpha_1=alpha, constrain_zero=False)['y_fit']
            BIC = get_stats_l1_model(y, fy, alpha)[-1]
            if BIC > BIC_REF:
                GO = False
            else:
                BIC_REF = BIC
                alpha_ref = alpha
        VERBOSE("alpha_ref: %f BIC_REF: %f" % (alpha_ref, BIC_REF))

        # search
        GO = True
        delta = 1.5
        alpha = alpha_ref
        while GO:
            alpha = alpha * delta
            fy = trend_filter(x, y, l_norm=1, alpha_1=alpha, constrain_zero=False)['y_fit']
            BIC = get_stats_l1_model(y, fy, alpha)[-1]
            if BIC > BIC_REF:
                GO = False
            else:
                BIC_REF = BIC
                alpha_ref = alpha
        VERBOSE("alpha_ref: %f BIC_REF: %f" % (alpha_ref, BIC_REF))
        GO = True
        delta = 1. / 1.5
        alpha = alpha_ref
        while GO:
            alpha = alpha * delta
            fy = trend_filter(x, y, l_norm=1, alpha_1=alpha, constrain_zero=False)['y_fit']
            BIC = get_stats_l1_model(y, fy, alpha)[-1]
            if BIC > BIC_REF:
                GO = False
            else:
                BIC_REF = BIC
                alpha_ref = alpha
        VERBOSE("alpha_ref: %f BIC_REF: %f" % (alpha_ref, BIC_REF))

        return trend_filter(x, y, l_norm=1, alpha_1=alpha_ref, constrain_zero=False)['y_fit']

    # new gts
    ngts = gts.copy()
    if isinstance(alpha,float):
        VERBOSE("Doing l1 trend filtering using alpha %f" % alpha)
        VERBOSE('North component')
        ngts.data[:, 1] = trend_filter(x, gts.data[:, 1] * 1.E3,l_norm=1, alpha_1=alpha, constrain_zero=False)['y_fit'] * 1.E-3
        VERBOSE('East  component')
        ngts.data[:, 2] = trend_filter(x, gts.data[:, 2] * 1.E3,l_norm=1, alpha_1=alpha, constrain_zero=False)['y_fit'] * 1.E-3
        VERBOSE('Up    component')
        ngts.data[:, 3] = trend_filter(x, gts.data[:, 3] * 1.E3,l_norm=1, alpha_1=alpha, constrain_zero=False)['y_fit'] * 1.E-3

    else:
        VERBOSE("Doing l1 trend filtering using optimal BIC criterion")
        VERBOSE('North component')
        ngts.data[:, 1] = f(x, gts.data[:, 1] * 1.E3) * 1.E-3
        VERBOSE('East  component')
        ngts.data[:, 2] = f(x, gts.data[:, 2] * 1.E3) * 1.E-3
        VERBOSE('Up    component')
        ngts.data[:, 3] = f(x, gts.data[:, 3] * 1.E3) * 1.E-3

    # return
    return ngts