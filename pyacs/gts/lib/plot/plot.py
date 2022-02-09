###############################################################################
def plot(self,
         title=None,
         loffset=True,
         loutliers=True,
         verbose=False,
         date=[],
         yaxis=None,
         min_yaxis=None,
         yupaxis=None,
         xticks_minor_locator=1,
         lcomponent=['N','E','U'],
         error_scale=1.0,
         lperiod=[[]],
         lvline=[],
         save_dir_plots='.',
         save=None,
         show=True,
         unit='mm',
         date_unit='cal',
         date_ref=0.0, center=True,
         superimposed=None,
         lcolor=['r','g','c','m','y','k','b'],
         label=None,
         legend=False,
         set_zero_at_date=None,
         grid=True,
         plot_size=None,
         info = [],
         xlabel_fmt = None,
         **kwargs):
###############################################################################
    """
    Create a plot of a North-East-Up time series and related info (offsets, outliers) using Matplotlib

    Coordinates of the time series are assumed to be in meters
    default plots units will be mm; Use unit='m' to get meters instead

    :param title: string to be added to the site name as a plot title
    :param loffset: boolean
        print a dash vertical line at offset dates
    :param loutliers: boolean
        print outliers
    :param verbose: boolean
        verbose mode
    :param date: [sdate,edate]
        start and end date for plots
        sdate and edate in decimal years if date_unit is either 'decyear' or 'cal', or in days if date_unit is 'days'
    :param yaxis: [min_y,max_y]
        min and max value for the yaxis
        if not provided automatically adjusted
    :param yupaxis: same as yaxis but applies to the up component only
    :param xticks_minor_locator: where xticks_minor_locator will be placed. Float when date_unit is 'decyear' or 'days', a string '%Y','%m','%d' is date_unit is 'cal'.  
    :param lcomponent: list of components to be plotted (default =['N','E','U'])
    :param error_scale: scaling factor for error bars (default = 1.0, 0 means no error bar)
    :param lperiod: list of periods to be drawn in background (color=light salmon)
    :param lvline: list of dates where vertical lines will be drawn in background (color=green)
    :param save_dir_plots: directory used for saving the plots
    :param save: name, save the plot into name, if simply True an automatic name is given
    :param show: boolean, is True show the plot
    :param unit: 'm','cm','mm', default='mm'
    :param date_unit: 'decyear' or 'cal' or 'days', default='decyear'
    :param date_ref: reference date, default=0.0
    :param center: boolean, if True the y_axis is centered around the mean value for the plotted period
    :param superimposed: if a gts is provided, it is superimposed to the master, default=None
    :param lcolor: color list used for the superimposed time series, default=['r','g','c','m','y','k','b']
    :param label: label for superimposed time series to be displayed in legend, default=None
    :param legend: boolean. Set true to display label for superimposed time series, default=False
    :param set_zero_at_date: date at which the master and superimposed gts will be equal (default=None). date can also be a list with [date,offset_north,offset_east,offset_up]
    :param plot_size: plot size as a tuple. Default, best guess.
    :param grid: boolean
    :param info: title to appear in time series subplots
    :param **kwargs: any argument to be passed to  matplotlib.pyplot.errorbar

    :note: The list of kwargs are:
    
    {  'linewidth' : 0,\
       'marker' : marker_main_symbol ,\
       'markersize' : marker_main_size, \
       'markerfacecolor' : marker_main_color,\
       'markeredgecolor' : marker_main_color,\
       'markeredgewidth' : marker_main_colorlw,\
       'ecolor' : error_bar_color,\
       'elinewidth' : error_bar_linewidth,\
       'capsize' : error_bar_capsize }

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


    ###########################################################################
    # lcolor
    ###########################################################################
    
    if ( superimposed is not None) and (isinstance(superimposed,list)) and (len(superimposed) >7): 
        lcolor = np.random.rand(len(superimposed),3)
    
    ###########################################################################
    # check data is not None
    ###########################################################################
    
    try:
        if self.data is None:
            # raise exception
            raise GtsInputDataNone(inspect.stack()[0][3],__name__,self)
    except GtsInputDataNone as error:
        # print PYACS WARNING
        print( error )
        return( self )


    ###########################################################################
    # ensure superimposed is either None or a list of gts
    ###########################################################################
    
    if superimposed is not None:
        if not isinstance(superimposed,list):
            superimposed=[superimposed] 

        for agts in superimposed:
            if not isinstance(agts,Gts):
                raise TypeError
    
    
    ###########################################################################
    # ensure lperiod is a list of list
    ###########################################################################
    lperiod = pyacs.lib.utils.__ensure_list_of_list( lperiod )
    
    ###########################################################################
    # init default plot settings
    ###########################################################################
    
    H_plot_settings , H_kwargs_errorbar , H_kwargs_errorbar_superimposed = pyacs.gts.lib.plot.init_plot_settings.__init_kwargs()

    ###########################################################################
    # init plot settings with user provided values
    ###########################################################################

    for k,v in kwargs.items():
        H_kwargs_errorbar[k] = v 
    
    # turn show to off by default
    plt.ioff()
    
    # some drawing default parameters
    params = {'legend.fontsize': 5}
    plt.rcParams.update(params)

    # For figure caption
    txt_component={}
    txt_component['N']='North'
    txt_component['E']='East'
    txt_component['U']='Up'

    # plot size
    if plot_size is None:
        if len(lcomponent)==1:plot_size = H_plot_settings['plot_size_1_component']
        if len(lcomponent)==2:plot_size = H_plot_settings['plot_size_2_component']
        if len(lcomponent)==3:plot_size = H_plot_settings['plot_size_3_component']


    lperiod_facecolor='r'
    lperiod_alpha=0.4


    ###########################################################################
    # format for printing duration in subplot titles
    ###########################################################################

    if date != []:duration=date[-1]-date[0]
    else:duration=self.data[-1,0]-self.data[0,0]
    
    ###########################################################################
    # x-tick label
    ###########################################################################
    # adjust tickers for the decyear case
    if date_unit=='decyear' or date_unit=='days':
        if xlabel_fmt != None:        
            fmt=ticker.FormatStrFormatter( xlabel_fmt )
        else:
            if duration >= 1.0:fmt=ticker.FormatStrFormatter('%4.0lf')
            if duration < 1.0:fmt=ticker.FormatStrFormatter('%6.2lf')
            if duration < 0.2:fmt=ticker.FormatStrFormatter('%6.3lf')
            if duration < 0.1:fmt=ticker.FormatStrFormatter('%6.4lf')
            
        plt.rcParams['axes.formatter.limits']=(-7,17)

    if date_unit=='cal':
        if xlabel_fmt is not None:        
            fmt = mdates.DateFormatter( xlabel_fmt )

    
    ###########################################################################
    # title and subtitle
    ###########################################################################
    
    if title is not None:
        stitle = title
    else:
        stitle = self.code

    ###########################################################################
    # info
    ###########################################################################
    
    if info == []:
        info = len(lcomponent) * ['']

    ###########################################################################
    # to_unit
    ###########################################################################

    if unit == 'cm':
        to_unit = 100.

    if unit == 'mm':
        to_unit = 1000.

    ###########################################################################
    # START PLOTTING
    ###########################################################################
    plt.ioff()
    f , ax = plt.subplots( len(lcomponent) , sharex=True , figsize=plot_size)


    ###########################################################################
    # rotate xticks labels
    ###########################################################################
    f.autofmt_xdate(bottom=0.2, rotation=30, ha='right')

# handle case where only one component and ax is not subscriptable
    if not isinstance(ax,np.ndarray):
        ax = np.array([ ax ])

    np_date , data = pyacs.gts.lib.plot.gts_to_ts.gts_to_ts( self , \
                                                             date=date, \
                                                             unit=unit , \
                                                             date_unit=date_unit , \
                                                             date_ref=date_ref, \
                                                             set_zero_at_date = set_zero_at_date, \
                                                             center = center )
   
# master time series
    ###########################################################################
    idx_ax = 0

    for component in lcomponent:
        # idx_component
        if component == 'N':
            idx = 0
        if component == 'E':
            idx = 1
        if component == 'U':
            idx = 2

        # subplot title
        str_title = pyacs.gts.lib.plot.make_stitle.make_stitle(component, info=info[ idx_ax ])
        ax[idx_ax].set_title(str_title, x=0, horizontalalignment='left', fontsize= '11')
   
        # plot the data
#        print( H_kwargs_errorbar )
        ax[idx_ax].errorbar(np_date,data[:,idx], yerr=data[:,idx+3]*error_scale, **H_kwargs_errorbar)

        # plot outliers
        if self.outliers != []:
            ax[idx_ax].plot(np_date[self.outliers],data[self.outliers,idx], 'ro' , markersize=H_kwargs_errorbar['markersize'])
            


        # Y-label    
        ax[idx_ax].set_ylabel(unit)

        # grid
        ax[ idx_ax ].grid( grid )

        # xticks
        if date_unit == 'decyear':
            ax[idx_ax].xaxis.set_major_formatter(fmt)

        if date_unit == 'cal' and xlabel_fmt is not None:
            ax[idx_ax].xaxis.set_major_formatter( fmt )

        idx_ax = idx_ax + 1

    # superimposed time series
    ###########################################################################
    
    if superimposed is None:
        wsuperimposed = []
    else:
        wsuperimposed = list( superimposed )
    
    # Get label
    
    if label is None:
        label = []
        for my_ts in wsuperimposed:
            label.append( my_ts.code ) 
    
    # test
    np_date = None
    data = None
    
    for i in np.arange( len(wsuperimposed) ):
        
        if verbose:
            print('-- plotting superimposed time series: ' , label[i])
    
        np_date , data = pyacs.gts.lib.plot.gts_to_ts.gts_to_ts( wsuperimposed[i] , \
                                                                 date=date, \
                                                                 unit=unit , \
                                                                 date_unit=date_unit , \
                                                                 date_ref=date_ref, \
                                                                 set_zero_at_date = set_zero_at_date, \
                                                                 center = center )
        

        idx_ax = 0
    
        for component in lcomponent:
            # idx_component
            if component == 'N':
                idx = 0
            if component == 'E':
                idx = 1
            if component == 'U':
                idx = 2
    
            # plot the data
            error_scale=0.
            ax[idx_ax].errorbar(np_date,data[:,idx], yerr=data[:,idx+3]*error_scale , color=lcolor[i] , label = label[i] )
            idx_ax = idx_ax +1 

    # plot offset_date of the master time series
    ###########################################################################
    if loffset:
        if date_unit == 'days':
            if (date_ref is None) or (date_ref == 0):
                # reference data is year.000
                date_ref=float(int(self.data[0,0]))
             
            ldate_offsets_in_date_unit = pyacs.lib.astrotime.decyear2mjd( self.offsets_dates ) - pyacs.lib.astrotime.decyear2mjd( date_ref ) 
    
        if date_unit == 'decyear':
            ldate_offsets_in_date_unit=np.array(self.offsets_dates)-date_ref
         
        if date_unit == 'cal':
            ldate_offsets_in_date_unit = pyacs.lib.astrotime.decyear2datetime(self.offsets_dates)
    
        idx_ax = 0
        for component in lcomponent:
            for odate in ldate_offsets_in_date_unit:
                ax[idx_ax].axvline(x= odate ,linewidth=2, color='k', linestyle='--')
            idx_ax = idx_ax + 1

    # lperiod option: periods to be highlighted 
    ###########################################################################
         
    if lperiod != [[]]:
        if date_unit == 'days':
            if (date_ref is None) or (date_ref == 0):
                # reference data is year.000
                date_ref=float(int(self.data[0,0]))

            lperiod=pyacs.lib.astrotime.day_since_decyear(np.array(lperiod).flatten(),date_ref).reshape(-1,2)
        
        if date_unit == 'decyear':
            lperiod=(np.array(lperiod).flatten()-date_ref).reshape(-1,2)
        
        if date_unit == 'cal':
            lperiod=np.array( pyacs.lib.astrotime.decyear2datetime(np.array(lperiod).flatten())).reshape(-1,2)
    
        # plot color background for periods
        idx_ax = 0
        for component in lcomponent:
            for period in lperiod:
                (sdate,edate)=period
                ax[idx_ax].axvspan(sdate, edate, facecolor=lperiod_facecolor, alpha=lperiod_alpha)
            
            idx_ax = idx_ax + 1
    
    # vertical lines option
    ###########################################################################

    if lvline != []:
    
        if date_unit == 'decyear':
            lvline_in_date_unit=np.array( lvline )-date_ref

        if date_unit == 'days':
            lvline_in_date_unit=np.array( lvline )-date_ref

        
        if date_unit == 'cal':
            lvline_in_date_unit = pyacs.lib.astrotime.decyear2datetime( lvline )
    
        idx_ax = 0
        for component in lcomponent:
            for odate in lvline_in_date_unit:
                ax[idx_ax].axvline(x= odate ,linewidth=2, color='g', linestyle='--')
            idx_ax = idx_ax + 1

        
    # title
    ###########################################################################
    plt.suptitle(stitle)

    # yaxis
    ###########################################################################


    # added JMN 13/01/2021 to set a minimum value for y-scale
    if (min_yaxis is not None) and (yaxis is None):

        idx_ax = 0
        for component in lcomponent:
            if component in ['N', 'E']:
                cyaxis=ax[idx_ax].get_ylim()
                ax[idx_ax].set_ylim( np.min([cyaxis[0],-min_yaxis]), np.max([cyaxis[1],min_yaxis]))
            if component == 'U':
                cyupaxis=ax[idx_ax].get_ylim()
                ax[idx_ax].set_ylim( np.min([cyupaxis[0],-min_yaxis]), np.max([cyupaxis[1],min_yaxis]))

            idx_ax = idx_ax + 1

    if yaxis is not None:
        idx_ax = 0
        for component in lcomponent:
            if component in ['N','E']:
                ax[idx_ax].set_ylim( yaxis[0] , yaxis[1] )
            if component == 'U' and yupaxis is not None:
                ax[idx_ax].set_ylim( yupaxis[0] , yupaxis[1] )
            
            idx_ax = idx_ax + 1


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
        idx_ax = 0
        for component in lcomponent:
            if component in ['N', 'E']:
                ax[idx_ax].set_xlim(np_date_x[0], np_date_x[1])
            if component == 'U' and yupaxis is not None:
                ax[idx_ax].set_xlim(np_date_x[0], np_date_x[1])

            idx_ax = idx_ax + 1

    # X Label
    ###########################################################################
    
    if date_unit=='decyear':
        if (date_ref is None) or (date_ref == 0.0):
            str_xlabel='year'
        else:
            str_xlabel=("decimal year since %s" %pyacs.lib.astrotime.decyear2datetime(date_ref))
    

    if date_unit=='days':
        
        if (date_ref is None) or (date_ref == 0):
            # reference data is year.000
            date_ref=float(int(self.data[0,0]))
        
        str_xlabel=("days since %s" %pyacs.lib.astrotime.decyear2datetime(date_ref))

    if date_unit=='cal':
        str_xlabel=("calendar date")

    
    ax[idx_ax -1 ].set_xlabel(str_xlabel)


    # xtcicks minor locator
    ###########################################################################    
    if date_unit == 'decyear':
        idx_ax = 0
        for component in lcomponent:
            ax[idx_ax].xaxis.set_minor_locator(MultipleLocator( float(xticks_minor_locator) ))
            idx_ax = idx_ax + 1

    if date_unit != 'days':
        idx_ax = 0
        if duration > 100:
            xticks_minor_locator = 10
        if duration > 1000:
            xticks_minor_locator = 100
        if duration > 10000:
            xticks_minor_locator = 1000
            
        for component in lcomponent:
            
            ax[idx_ax].xaxis.set_minor_locator(MultipleLocator( float(xticks_minor_locator) ))
            idx_ax = idx_ax + 1


    if date_unit == 'cal':
        # default
        if self.data[-1,0] - self.data[0,0] > 5.:
            locator = mdates.YearLocator()
        else:
            locator = mdates.MonthLocator()
        
        # user provided default    
        if xticks_minor_locator == '%Y':
            locator = mdates.YearLocator()   # every year
        if xticks_minor_locator == '%m':
            locator = mdates.MonthLocator()   # every month
        if xticks_minor_locator == '%d':
            locator = mdates.DayLocator()   # every day
        
        idx_ax = 0
        for component in lcomponent:
            ax[idx_ax].xaxis.set_minor_locator( locator )
            idx_ax = idx_ax + 1
        

    # legend
    ###########################################################################

    if legend:
        idx_ax = 0
        for component in lcomponent:
            ax[idx_ax].legend()
            idx_ax = idx_ax + 1
        

    # tight_layout
    ###########################################################################

    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace= H_plot_settings['hspace'])
    
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
                fname = os.path.normpath(save_dir_plots+'/'+save)

        elif save:
            # save is simply set to True for automatic output file naming
                fname = os.path.normpath( save_dir_plots + '/' + self.code+'.png' )
        
        # saving output file
        
        if verbose:
            print("-- Saving file to ",fname)
        
        # creates the directory of it does not exist
        
        if os.path.sep in fname:
            os.makedirs(os.path.dirname(fname), exist_ok=True)
        
#        plt.savefig(fname, dpi=150, facecolor='w', edgecolor='w', orientation='portrait', papertype='a4', format='png',transparent=False, pad_inches=0.1)
# change 05/04/2021 papertype appears to be decommissioned for matplotlib 3.3
        plt.savefig(fname, dpi=150, facecolor='w', edgecolor='w', orientation='portrait', format='png',transparent=False, pad_inches=0.1)

    # show
    ###########################################################################
    if show: 
        # do not block program execution when displaying figures
        plt.ion()
        plt.show()
    else:
        plt.close( f )

    # returns the Gts for pyacs convention compatibility        
    new_Gts=self.copy()
    
    return( new_Gts )
    
