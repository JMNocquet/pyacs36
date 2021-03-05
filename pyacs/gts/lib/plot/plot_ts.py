###############################################################################
def plot_ts(ax, np_date, data, yerr,
            loutliers=[],
            yaxis=None,
            error_scale=1.0,
            **H_kwargs):
    ###############################################################################

    # import
    import pyacs.lib.astrotime
    import numpy as np

    if yaxis is not None:
        (ymin, ymax) = yaxis
        ax.set_ylim(ymin, ymax)

    # plot data
    ax.errorbar(np_date, data, yerr=yerr * error_scale, **H_kwargs)
    ax.grid(True)


    return
