###############################################################################
def plot_component(self,
         component='E',
         lorder=[],
         shift = 10.,
         Hshift={},
         title=None,
         loffset=True,
         loutliers=True,
         verbose=False,
         date=[],
         yaxis=None,
         min_yaxis=None,
         yupaxis=None,
         xticks_minor_locator=1,
         error_scale=1.0,
         lperiod=[[]],
         lvline=[],
         save_dir_plots='.',
         save=None,
         show=True,
         unit='mm',
         date_unit='cal',
         date_ref=0.0,
         set_zero_at_date=None,
         superimposed=None,
         lcolor=[],
         label=None,
         legend=0.2,
         grid=True,
         plot_size=None,
         xlabel_fmt=None,
         **kwargs):
    ###############################################################################
    """Plot one component for multiple time series using Matplotlib.

    Coordinates in meters; default display unit mm (unit='m' for meters).

    Parameters
    ----------
    component : str, optional
        Component to plot ('N', 'E', 'U'). Default is 'E'.
    lorder : list, optional
        Site code order (top to bottom). Default is [].
    shift : float, optional
        Vertical shift between series in mm. Default is 10.
    Hshift : dict, optional
        Per-code offset {'code': offset}. Default is {}.
    title : str, optional
        Plot title. Default is None.
    loffset : bool, optional
        Draw offset lines. Default is True.
    loutliers : bool, optional
        Plot outliers. Default is True.
    verbose : bool, optional
        Verbose mode. Default is False.
    date : list, optional
        [sdate, edate] for range. Default is [].
    yaxis, min_yaxis, yupaxis : list, optional
        Y-axis limits; auto if not set.
    xticks_minor_locator : float or str, optional
        Minor tick spacing. Default is 1.
    error_scale : float, optional
        Error bar scale (0 = none). Default is 1.0.
    lperiod : list, optional
        Background periods. Default is [[]].
    lvline : list, optional
        Dates for vertical lines. Default is [].
    save_dir_plots : str, optional
        Save directory. Default is '.'.
    save : bool or str, optional
        Save figure. Default is None.
    show : bool, optional
        Show plot. Default is True.
    unit : str, optional
        'm', 'cm', 'mm'. Default is 'mm'.
    date_unit : str, optional
        'decyear', 'cal', 'days'. Default is 'cal'.
    date_ref : float, optional
        Reference date. Default is 0.0.
    set_zero_at_date : float, optional
        Date to align series. Default is None.
    superimposed : Gts, optional
        Overlay Gts. Default is None.
    lcolor : list, optional
        Colors for superimposed. Default is [].
    label : str, optional
        Legend label. Default is None.
    legend : float, optional
        Legend space (e.g. 0.2 = 20%). Default is 0.2.
    grid : bool, optional
        Draw grid. Default is True.
    plot_size : tuple, optional
        Figure size. Default is None.
    xlabel_fmt : str, optional
        X-axis label format. Default is None.
    **kwargs : dict, optional
        Passed to matplotlib.pyplot.errorbar.

    Returns
    -------
    None
    """

    ###########################################################################
    # import
    ###########################################################################

    # system
    import inspect
    import numpy as np
    import os.path

    # matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    from matplotlib.ticker import MultipleLocator
    import matplotlib.dates as mdates

    # pyacs
    from pyacs.gts.lib.errors import GtsInputDataNone
    import pyacs.gts.lib.plot.init_plot_settings
    import pyacs.lib.astrotime
    from pyacs.gts.Gts import Gts
    import pyacs.lib.utils
    import pyacs.gts.lib.plot.make_stitle
    import pyacs.gts.lib.plot.gts_to_ts

    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG

    import inspect
    VERBOSE("Running Sgts.%s" % inspect.currentframe().f_code.co_name)

    # debug
    from icecream import ic

    # Validate input parameters
    if component not in ['N', 'E', 'U']:
        raise ValueError(f"Invalid component: {component}. Must be one of ['N', 'E', 'U']")

    if unit not in ['m', 'cm', 'mm']:
        raise ValueError(f"Invalid unit: {unit}. Must be one of ['m', 'cm', 'mm']")

    if date_unit not in ['decyear', 'cal', 'days']:
        raise ValueError(f"Invalid date_unit: {date_unit}. Must be one of ['decyear', 'cal', 'days']")

    if not isinstance(error_scale, (int, float)) or error_scale < 0:
        raise ValueError(f"Invalid error_scale: {error_scale}. Must be a non-negative number")

    if not isinstance(shift, (int, float)) or shift <= 0:
        raise ValueError(f"Invalid shift: {shift}. Must be a positive number")

    if not isinstance(legend, (int, float)) or not (0 < legend < 1):
        raise ValueError(f"Invalid legend: {legend}. Must be between 0 and 1")

    if date and len(date) != 2:
        raise ValueError("date must be a list of two elements [start_date, end_date]")

    if date_ref is not None and not isinstance(date_ref, (int, float)):
        raise ValueError("date_ref must be a number")

    # lcolor
    if lcolor == []:
        import matplotlib.colors as mcolors
        lcolor = list(mcolors.BASE_COLORS.keys())[:-2] + list(mcolors.TABLEAU_COLORS.keys()) + list(mcolors.XKCD_COLORS.keys())

    ###########################################################################
    # ensure lperiod is a list of list
    ###########################################################################
    lperiod = pyacs.lib.utils.__ensure_list_of_list(lperiod)

    ###########################################################################
    # init default plot settings
    ###########################################################################

    H_plot_settings, H_kwargs_errorbar, H_kwargs_errorbar_superimposed = pyacs.gts.lib.plot.init_plot_settings.__init_kwargs()

    ###########################################################################
    # init plot settings with user provided values
    ###########################################################################

    for k, v in kwargs.items():
        H_kwargs_errorbar[k] = v

    # turn show to off by default
    plt.ioff()

    # some drawing default parameters
    params = {'legend.fontsize': 11}
    plt.rcParams.update(params)

    # For figure caption
    txt_component = {}
    txt_component['N'] = 'North'
    txt_component['E'] = 'East'
    txt_component['U'] = 'Up'

    # plot size
    if plot_size is None:
        plot_size = H_plot_settings['plot_size_3_component']

    lperiod_facecolor = 'r'
    lperiod_alpha = 0.4

    ###########################################################################
    # title and subtitle
    ###########################################################################

    if title is not None:
        stitle = title
    else:
        stitle = txt_component[component]

    ###########################################################################
    # to_unit
    ###########################################################################

    if unit == 'cm':
        to_unit = 100.
    elif unit == 'mm':
        to_unit = 1000.
    else:
        to_unit = 1.

    ###########################################################################
    # START PLOTTING
    ###########################################################################
    plt.ioff()

    ###########################################################################
    # handle order
    ###########################################################################

    if lorder == []:
        # use alpha_numeric order
        lorder = sorted(self.lcode())

    lorder.reverse()

    ###########################################################################
    # format for printing duration in subplot titles
    ###########################################################################

    if date != []:
        duration = date[-1] - date[0]
    else:
        try:
            duration = self.__dict__[lorder[0]].data[-1, 0] - self.__dict__[lorder[0]].data[0, 0]
        except (KeyError, IndexError, AttributeError) as e:
            raise ValueError(f"Could not determine duration from time series: {str(e)}")

    ###########################################################################
    # x-tick label
    ###########################################################################
    # adjust tickers for the decyear case
    if date_unit == 'decyear' or date_unit == 'days':
        if xlabel_fmt != None:
            fmt = ticker.FormatStrFormatter(xlabel_fmt)
        else:
            if duration >= 1.0: fmt = ticker.FormatStrFormatter('%4.0lf')
            if duration < 1.0: fmt = ticker.FormatStrFormatter('%6.2lf')
            if duration < 0.2: fmt = ticker.FormatStrFormatter('%6.3lf')
            if duration < 0.1: fmt = ticker.FormatStrFormatter('%6.4lf')

        plt.rcParams['axes.formatter.limits'] = (-7, 17)

    if date_unit == 'cal':
        if xlabel_fmt is not None:
            fmt = mdates.DateFormatter(xlabel_fmt)

    ###########################################################################
    # handle shift
    ###########################################################################

    for i,code in enumerate(lorder):
        if code in Hshift.keys():
            Hshift[code] = shift * i + Hshift[code]
        else:
            Hshift[code] = shift * i

    ###########################################################################
    # init figure
    ###########################################################################
    Hctoidx={'N':0,'E':1,'U':2}
    idxc = Hctoidx[component]

    f, ax = plt.subplots(1,figsize=plot_size)

    # plot settings
    # grid
    ax.grid(grid)
    # xticks
    if date_unit == 'decyear':
        ax.xaxis.set_major_formatter(fmt)

    if date_unit == 'cal' and xlabel_fmt is not None:
        ax.xaxis.set_major_formatter(fmt)

    # subplot title
    ax.set_title(stitle, x=0, horizontalalignment='left', fontsize='11')

    ###########################################################################
    # rotate xticks labels
    ###########################################################################
    f.autofmt_xdate(bottom=0.2, rotation=30, ha='right')

    ###########################################################################
    # date & unit transformation
    ###########################################################################

    ts = self.sub(linclude=lorder)
    if date != []:
        ts = ts.gts('extract_periods',date)

    if set_zero_at_date is None:
        min_date,max_date=ts.get_dates(fmt='decyear')
        set_zero_at_date = (min_date + max_date) / 2

    yaxismax = -99
    yaxismin = 99

    for i,code in enumerate(lorder):
        try:
            # Validate Gts object
            if code not in self.__dict__:
                WARNING(f"Time series {code} not found in Sgts object")
                continue

            gts = self.__dict__[code]


            if not isinstance(gts, Gts):
                WARNING(f"Invalid Gts object for {code}")
                continue

            if gts.data is None or len(gts.data) == 0:
                WARNING(f"No data available for {code}")
                continue

            # Convert Gts to time series
            try:
                np_date, data = pyacs.gts.lib.plot.gts_to_ts.gts_to_ts(
                    gts,
                    date=date,
                    unit=unit,
                    date_unit=date_unit,
                    set_zero_at_date=set_zero_at_date,
                    center=False
                )
            except Exception as e:
                WARNING(f"Failed to convert Gts {code} to time series: {str(e)}")
                continue

            # Validate converted data
            if np_date is None or data is None or len(np_date) == 0 or len(data) == 0:
                WARNING(f"No valid data after conversion for {code}")
                continue

            # plot the data
            H_kwargs_errorbar['markerfacecolor'] = lcolor[i]
            ax.errorbar(np_date, data[:, idxc] + Hshift[code], 
                       yerr=data[:, idxc + 3] * error_scale, 
                       label=code, 
                       **H_kwargs_errorbar)
#            ax.errorbar(np_date, data[:, idxc] - np.median(data[:, idxc]) + Hshift[code], 
#                       yerr=data[:, idxc + 3] * error_scale, 
#                       label=code, 
#                       **H_kwargs_errorbar)

            # superimposed time series
            if superimposed is not None and code in superimposed.lcode():
                try:
                    sup_gts = superimposed.__dict__[code]
                    if not isinstance(sup_gts, Gts):
                        WARNING(f"Invalid superimposed Gts object for {code}")
                        continue

                    sup_np_date, sup_data = pyacs.gts.lib.plot.gts_to_ts.gts_to_ts(
                        sup_gts,
                        date=date,
                        unit=unit,
                        date_unit=date_unit,
                        set_zero_at_date=set_zero_at_date,
                        center=False
                    )

                    if sup_np_date is not None and sup_data is not None:
                        ax.errorbar(sup_np_date, 
                                  sup_data[:, idxc] - np.median(sup_data[:, idxc]) + Hshift[code], 
                                  yerr=sup_data[:, idxc + 3] * error_scale, 
                                  color=lcolor[i])
                except Exception as e:
                    WARNING(f"Failed to plot superimposed time series for {code}: {str(e)}")

            if loffset:
                try:
                    if date_unit == 'days':
                        if (date_ref is None) or (date_ref == 0):
                            date_ref = float(int(gts.data[0, 0]))
                        ldate_offsets_in_date_unit = pyacs.lib.astrotime.decyear2mjd(
                            gts.offsets_dates) - pyacs.lib.astrotime.decyear2mjd(date_ref)
                    elif date_unit == 'decyear':
                        ldate_offsets_in_date_unit = np.array(gts.offsets_dates) - date_ref
                    else:  # date_unit == 'cal'
                        ldate_offsets_in_date_unit = pyacs.lib.astrotime.decyear2datetime(gts.offsets_dates)

                    for odate in ldate_offsets_in_date_unit:
                        ax.axvline(x=odate, linewidth=2, color='k', linestyle='--')
                except Exception as e:
                    WARNING(f"Failed to plot offsets for {code}: {str(e)}")

        except Exception as e:
            WARNING(f"Error processing time series {code}: {str(e)}")
            continue

    # lperiod option: periods to be highlighted
    ###########################################################################

    if lperiod != [[]]:
        if date_unit == 'days':
            if (date_ref is None) or (date_ref == 0):
                # reference data is year.000
                date_ref = float(int(self.data[0, 0]))

            lperiod = pyacs.lib.astrotime.day_since_decyear(np.array(lperiod).flatten(), date_ref).reshape(-1, 2)

        if date_unit == 'decyear':
            lperiod = (np.array(lperiod).flatten() - date_ref).reshape(-1, 2)

        if date_unit == 'cal':
            lperiod = np.array(pyacs.lib.astrotime.decyear2datetime(np.array(lperiod).flatten())).reshape(-1, 2)

        # plot color background for periods
        for period in lperiod:
            (sdate, edate) = period
            ax.axvspan(sdate, edate, facecolor=lperiod_facecolor, alpha=lperiod_alpha)

    # vertical lines option
    ###########################################################################

    if lvline != []:

        if date_unit == 'decyear':
            lvline_in_date_unit = np.array(lvline) - date_ref

        if date_unit == 'days':
            lvline_in_date_unit = np.array(lvline) - date_ref

        if date_unit == 'cal':
            lvline_in_date_unit = pyacs.lib.astrotime.decyear2datetime(lvline)

        for odate in lvline_in_date_unit:
            ax.axvline(x=odate, linewidth=2, color='g', linestyle='--')

    # title
    ###########################################################################
    # plt.suptitle(stitle)

    # date
    ###########################################################################

    if date != []:

        # case 'days'
        if date_unit == 'days':
            if (date_ref is None) or (date_ref == 0):
                np_date_x = pyacs.lib.astrotime.decyear2mjd(np.array(date)) - pyacs.lib.astrotime.decyear2mjd(
                    date[0])
            else:
                np_date_x = pyacs.lib.astrotime.decyear2mjd(np.array(date)) - pyacs.lib.astrotime.decyear2mjd(
                    date_ref)

        # case 'decyear'
        if date_unit == 'decyear':
            if (date_ref is None) or (date_ref == 0):
                np_date_x = np.array(date)
            else:
                np_date_x = np.array(date) - date_ref

        # case 'cal'
        if date_unit == 'cal':
            np_date_x = pyacs.lib.astrotime.decyear2datetime(np.array(date))

        ax.set_xlim(np_date_x[0], np_date_x[1])

    # X Label
    ###########################################################################

    if date_unit == 'decyear':
        if (date_ref is None) or (date_ref == 0.0):
            str_xlabel = 'year'
        else:
            str_xlabel = ("decimal year since %s" % pyacs.lib.astrotime.decyear2datetime(date_ref))

    if date_unit == 'days':

        if (date_ref is None) or (date_ref == 0):
            # reference data is year.000
            date_ref = float(int(self.data[0, 0]))

        str_xlabel = ("days since %s" % pyacs.lib.astrotime.decyear2datetime(date_ref))

    if date_unit == 'cal':
        str_xlabel = ("calendar date")

    ax.set_xlabel(str_xlabel)

    # xtcicks minor locator
    ###########################################################################
    if date_unit == 'decyear':
        ax.xaxis.set_minor_locator(MultipleLocator(float(xticks_minor_locator)))

    if date_unit != 'days':
        idx_ax = 0
        if duration > 100:
            xticks_minor_locator = 10
        if duration > 1000:
            xticks_minor_locator = 100
        if duration > 10000:
            xticks_minor_locator = 1000

        ax.xaxis.set_minor_locator(MultipleLocator(float(xticks_minor_locator)))

    if date_unit == 'cal':
        # default
        if self.__dict__[code].data[-1, 0] - self.__dict__[code].data[0, 0] > 5.:
            locator = mdates.YearLocator()
        else:
            locator = mdates.MonthLocator()

        # user provided default
        if xticks_minor_locator == '%Y':
            locator = mdates.YearLocator()  # every year
        if xticks_minor_locator == '%m':
            locator = mdates.MonthLocator()  # every month
        if xticks_minor_locator == '%d':
            locator = mdates.DayLocator()  # every day

        ax.xaxis.set_minor_locator(locator)

    # y-axis label

    ax.set_ylabel(unit)


    # legend
    ###########################################################################

    handles, labels = ax.get_legend_handles_labels()

    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * (1-legend), box.height])

    # Put a legend to the right of the current axis
    ax.legend(handles[::-1], labels[::-1],loc='center left', bbox_to_anchor=(1, .5),markerscale=3. )

    # adjust yxais
    ###########################################################################

    tmpts = self.gts('extract_periods',date)
    tmpts = tmpts.gts('set_zero_at_date',set_zero_at_date)

    yaxismax=-99999
    yaxismin=99999

    for code in tmpts.lcode():
        try:
            maxx = np.max(tmpts.__dict__[code].data[:, idxc+1])*1E3 + Hshift[code] + shift  # pyright: ignore[reportUndefinedVariable]
            minx = np.min(tmpts.__dict__[code].data[:, idxc+1])*1E3 + Hshift[code] - shift
            #print(code, maxx, minx)
            yaxismax = np.max([yaxismax, maxx])
            yaxismin = np.min([yaxismin, minx])
        except:
            continue

    ax.set_ylim(yaxismin, yaxismax)

    # tight_layout
    ###########################################################################

    #plt.tight_layout()

    #plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=H_plot_settings['hspace'])

    # end of subplots
    ###########################################################################

    # save file as png if save option is set
    ###########################################################################

    if save is not None:

        # case save is a string

        if isinstance(save, str):

            # check whether save includes a path or simply a file name

            if os.path.sep in save:
                # save includes a path information
                fname = save
            else:
                # else
                fname = os.path.normpath(save_dir_plots + '/' + save)

        elif save:
            # save is simply set to True for automatic output file naming
            fname = os.path.normpath(save_dir_plots + '/' + self.code + '.png')

        # saving output file

        if verbose:
            print("-- Saving file to ", fname)

        # creates the directory of it does not exist

        if os.path.sep in fname:
            os.makedirs(os.path.dirname(fname), exist_ok=True)

        #        plt.savefig(fname, dpi=150, facecolor='w', edgecolor='w', orientation='portrait', papertype='a4', format='png',transparent=False, pad_inches=0.1)
        # change 05/04/2021 papertype appears to be decommissioned for matplotlib 3.3
        plt.savefig(fname, dpi=150, facecolor='w', edgecolor='w', orientation='portrait', format='png',
                    transparent=False, pad_inches=0.1)

    # show
    ###########################################################################
    if show:
        # do not block program execution when displaying figures
        plt.ion()
        plt.show()
    else:
        plt.close(f)

    # returns the Sgts for pyacs convention compatibility
    return ( self )

