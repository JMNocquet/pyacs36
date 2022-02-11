def plot_data_sum(self, period=None):
    """
    generate an interactive plot showing available/missing data

    param period: period used to checkk data availability. Uses pandas syntax. See examples.
    return self
    examples:
            ts.plot_data_sum(period='2008') willl show a summary for the whole year 2008
            ts.plot_data_sum(period=['2008-01-15','2011-01-01']) willl show a summary for the period
            Dispay will include all time series in ts. Combine with ts selection to obtain more specific behaviour.
            ts.sub(linclude=['XXXX','YYYY']).plot_data_sum(period='2008') will plot data availability only for 'XXXX' and 'YYYY'
            ts.sel_radius('XXXX',[0,50]).sel_period([2008,2015],min_data=1000).plot_data_sum(period='2008') will show data availability for year 2008, only for sites less than 50 km from XXXX and having a minimum of 1000 days between 2008 and 2015.
    """
    import pandas
    import matplotlib.pyplot as plt
    import pyacs.message.verbose_message as VERBOSE
    import numpy as np
    import pyacs.lib.astrotime as at

    VERBOSE("Gts to tensor conversion for %d sites" % self.n())
    T, np_names, np_seconds_sol = self.to_obs_tensor(rounding='day')
    df = pandas.DataFrame(T[:, :, 0])
    df.index = at.seconds2datetime(np_seconds_sol)

    # decipher period
    if isinstance(period, list):
        a = df.loc[period[0]:period[1]]
    else:
        a = df.loc[period]

    npa = np.copy(a.to_numpy())
    # revert nan
    lidx_nan = np.where(np.isnan(npa))
    lidx_nnan = np.where(~np.isnan(npa))
    npa[lidx_nan] = 0
    npa[lidx_nnan] = np.nan
    npa = npa + np.arange(npa.shape[1] - 1, -1, -1)
    bnpa = np.ones(npa.shape) * np.arange(npa.shape[1])

    VERBOSE("Making plot")

    fig, ax = plt.subplots(figsize=(7, 9))

    ax.plot(a.index, bnpa, color='lime', linewidth=3)
    ax.plot(a.index, npa, color='red', linewidth=3)
    ax.set_yticks(np.arange(a.shape[1]))
    ax.set_yticklabels(np.flip(np_names))
    ax.tick_params(axis='y', which='major', labelsize=6)
    # ax.set_xlim(a.index[0],a.index[1])
    ax.set_ylim(-0.5, npa.shape[1] - 0.5)
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45)

    return self
