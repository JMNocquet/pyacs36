###############################################################################
## Gts Plot method using matplotlib
###############################################################################

import numpy as np
from pyacs.gts.Gts import Gts
import inspect


# printing options for debug
np.set_printoptions(precision=3,threshold=10000,linewidth=150)


###############################################################################
# GLOBAL PARAMETERS FOR PLOTS         
###############################################################################

# FIGURE SIZE

plot_size_1_component=(10,3.5)
plot_size_2_component=(8,8)
plot_size_3_component=(7,8)

# SPACE BETWEEN SUBPLOTS
hspace=0.3

# PLOT SETTINGS FOR MAIN MAKERS,ERRORBARS AND SUPERIMPOSED  
marker_main_symbol='o'
marker_main_color='b'
marker_main_size=4
marker_main_colorlw=0.0
marker_main_linewidth=0
marker_main_linecolor='k'

error_bar_color='grey'
error_bar_linewidth=1
error_bar_capsize=0

line_2_color='r'
line_2_width=.5
line_2_style='-'
marker_2_symbol=None
marker_2_color='r'
marker_2_size=None
marker_2_colorlw=0.0
marker_2_linewidth=1
marker_2_linecolor='r'

H_kwargs_plot_date={}
H_kwargs_plot_date_superimposed={}
H_kwargs_errorbar={}
H_kwargs_errorbar_superimposed={}


def __init_kwargs():

    global H_kwargs_plot_date
    global H_kwargs_plot_date_superimposed
    global H_kwargs_errorbar
    global H_kwargs_errorbar_superimposed

    
    H_kwargs_plot_date={\
                        'marker' : marker_main_symbol,\
                        'markersize' : marker_main_size, \
                        'markerfacecolor' : marker_main_color,\
                        'markeredgecolor' : marker_main_color,\
                        'markeredgewidth' : marker_main_colorlw }
    
    H_kwargs_plot_date_superimposed={ \
                        'color' : line_2_color,\
                        'linewidth' : line_2_width,\
                        'linestyle' : line_2_style,\
                        'marker' : marker_2_symbol,\
                        'markersize' : marker_2_size, \
                        'markerfacecolor' : marker_2_color,\
                        'markeredgecolor' : marker_2_color,\
                        'markeredgewidth' : marker_2_colorlw }
    
    
    H_kwargs_errorbar={'linewidth' : 0,\
                       'marker' : marker_main_symbol ,\
                       'markersize' : marker_main_size, \
                       'markerfacecolor' : marker_main_color,\
                       'markeredgecolor' : marker_main_color,\
                       'markeredgewidth' : marker_main_colorlw,\
                       'ecolor' : error_bar_color,\
                       'elinewidth' : error_bar_linewidth,\
                       'capsize' : error_bar_capsize }
    
    H_kwargs_errorbar_superimposed={ \
                        'color' : line_2_color,\
                        'linewidth' : line_2_width,\
                        'linestyle' : line_2_style,\
                        'marker' : marker_2_symbol,\
                        'markersize' : marker_2_size, \
                        'markerfacecolor' : marker_2_color,\
                        'markeredgecolor' : marker_2_color,\
                        'markeredgewidth' : marker_2_colorlw, \
                        'ecolor' : error_bar_color,\
                        'elinewidth' : error_bar_linewidth,\
                        'capsize' : error_bar_capsize }

# OTHER PLOT SETTINGS
# lperiod color
lperiod_facecolor='r'
lperiod_alpha=0.4


###############################################################################
# SOME ROUTINES         
###############################################################################

###############################################################################
def __ensure_list_of_list(ll):
###############################################################################
    """
    Ensures ll is a list of lists
    e.g. [a,b] returns [[a,b]],[[a,b]] returns [[a,b]]
    """

    if not isinstance(ll[0],list):
        return([ll])
    else:
        return(ll)

###############################################################################
def __window_ts(data,xmin=None,xmax=None):
###############################################################################
    if (xmin is None):xmin=-1
    if (xmax is None):xmax=1E20
    
    lindex=np.where( (data[:,0] > xmin) & (data[:,0] < xmax) )

    return(data[lindex])
    

###############################################################################
# DATES ROUTINES         
###############################################################################

###############################################################################
def __decyear2num(dates):
###############################################################################
    """
    converts decimal year to numerical dates used by matplotlib.dates
    """
    
    import pyacs.lib.astrotime

    import matplotlib.dates
    import datetime
    new_dates=[]
    for date in dates:
        year=int(date)
        (mday,month,ut)=pyacs.lib.astrotime.decyear2cal(date)
        t=datetime.date(year,month, mday)
        new_dates.append(matplotlib.dates.date2num(t)+ut)
    return(new_dates)


###############################################################################
# OTHER SUBROUTINES
###############################################################################

###############################################################################
def __make_stitle(component,velocity,data,unit, info=[]):
###############################################################################
    """
    
    :param info: a list of keywords of information to be displayed: 'wrms','bias','vel','period'
    """
    
    
    import pyacs.lib.astrotime

    str_title=''
  
    # component
    if component=='N':txt_component='North'
    if component=='E':txt_component='East'
    if component=='U':txt_component='Up'
    
    str_title = str_title+("%s " % (txt_component) )
    
    # velocity
    if 'vel' in info:
        if isinstance(velocity,float):
            str_title = str_title+(" V = %.2lf %s/yr" % (velocity, unit) )

    # wrms
    if 'wrms' in info:
        wmean= np.sum( (data[:,1]/data[:,2]**2 ) ) / np.sum(1./data[:,2]**2)
        wrms=np.sqrt( np.sum( ( (data[:,1]-wmean)**2/data[:,2]**2 ) ) / np.sum(1./data[:,2]**2) )
        str_title = str_title+(" wrms=%5.1lf " % (wrms) )
    
    # period  
    if 'period' in info:
        sdoy,_ut=pyacs.lib.astrotime.decyear2dayno(data[0,0])
        syear=int(data[0,0])
        
        edoy,_ut=pyacs.lib.astrotime.decyear2dayno(data[-1,0])
        eyear=int(data[-1,0])
        
        duration=data[-1,0]-data[0,0]
        
        str_period=(" %04d/%03d - %04d/%03d - %4.1lf yr" % (syear,sdoy,eyear,edoy,duration))
        str_title=str_title + (" #%d (%s)" % ( data.shape[0], str_period )) 

    # bias
    if 'bias' in info:
        bias=np.median((data[:,1]))
        str_title=str_title+" bias=%5.1lf " % (bias)

    return(str_title)

###############################################################################
# PLOT SUBROUTINES
###############################################################################

###############################################################################
def __plot_lperiod(plt,lperiod,date_unit,ref_date):
###############################################################################

    lperiod = __ensure_list_of_list(lperiod)

    if lperiod == [[]]:return

    import pyacs.lib

    if date_unit == 'days':
        lperiod=pyacs.lib.astrotime.day_since_decyear(np.array(lperiod).flatten(),ref_date).reshape(-1,2)
    
    if date_unit == 'decyear':
        
        lperiod=(np.array(lperiod).flatten()-ref_date).reshape(-1,2)
    
    if date_unit == 'cal':
        lperiod=np.array( __decyear2num(np.array(lperiod).flatten())).reshape(-1,2)

    # plot color background for periods
    for period in lperiod:
        (sdate,edate)=period
        plt.axvspan(sdate, edate, facecolor=lperiod_facecolor, alpha=lperiod_alpha)
    
    return 

###############################################################################
def __plot_lvline(plt,lvline,date_unit,ref_date):
###############################################################################
    
    import pyacs.lib
    
    if date_unit == 'days':
        ldate=pyacs.lib.astrotime.day_since_decyear(lvline,ref_date)
    if date_unit == 'decyear':
        ldate=lvline
    if date_unit == 'cal':
        ldate = __decyear2num(lvline)


    for odate in ldate:
            plt.axvline(x= odate ,linewidth=2, color='g', linestyle='--')

    
    for date in lvline:
        plt.axvline(x=float(date),linewidth=2, color='r', linestyle='--')



###############################################################################
def __plot_ts( plt, data, \
            loutliers=[], \
            date=[], date_unit='decyear', ref_date=None, \
            yaxis=None, \
            center=True, \
            superimposed=False, \
            verbose=False, \
            error_scale=1.0,\
            **H_kwargs):
###############################################################################

    # import
    import pyacs.lib.astrotime


    # handle ref_date
    
    if date_unit=='decyear':
        if ref_date is None:
            # reference data is 0.0
            ref_date=0.0
        
        data[:,0]=data[:,0]-ref_date
        
    if date_unit=='days':
        
        if (ref_date is None) or (ref_date == 0):
            # reference data is year.000
            ref_date=float(int(data[0,0]))
        
        dates=data[:,0]
        new_dates=pyacs.lib.astrotime.day_since_decyear(dates,ref_date)
        data[:,0]=np.array(new_dates)

    
    # Get xmax and xmin 
    
    if (date is  not None) and (date != []):
        
#        if date_unit=='cal':
#            xmin = __decyear2num([date[0]])
#            xmax = __decyear2num([date[1]])
        
        
        if date_unit=='decyear':
            xmin=date[0]-ref_date;xmax=date[1]-ref_date

        if date_unit=='days':
            xmin = pyacs.lib.astrotime.day_since_decyear(date[0] , ref_date )
            xmax = pyacs.lib.astrotime.day_since_decyear(date[1] , ref_date )
            
    else:

        if date_unit=='cal':
#            xmin = data[0 ,0] - 0.5
#            xmax = data[-1,0] + 0.5
            xmin = data[0 ,0] 
            xmax = data[-1,0] 

        
        if date_unit=='decyear':
            one_day_dec_year=1/365.25 
            xmin=data[0,0]-one_day_dec_year
            xmax=data[-1,0]+one_day_dec_year

        if date_unit=='days':
            one_day=1. 
            xmin=data[0,0]-one_day
            xmax=data[-1,0]+one_day
    
    
    # slice and center ts

    centered_data = __window_ts(data,xmin=xmin,xmax=xmax)

    if not isinstance(centered_data, np.ndarray):
        print('! WARNING: problems with data. Will not create the plot')
        return()
    
    
    if center:
        
        yaxis_shift = np.median(centered_data[:,1])

        centered_data[:,1]=centered_data[:,1]-yaxis_shift
        
    else:

        yaxis_shift = 0
        
    
    if date_unit in ['decyear','days']:
        plt.xlim(xmin=xmin,xmax=xmax)
    
    if yaxis:
        (ymin,ymax)=yaxis
        plt.ylim(ymin=ymin,ymax=ymax)
    

    # plot data

    if date_unit=='cal':
        plt.errorbar(pyacs.lib.astrotime.decyear2datetime(centered_data[:,0]),centered_data[:,1], yerr=centered_data[:,2]*error_scale, **H_kwargs)
            
    if date_unit=='decyear':
        plt.errorbar(centered_data[:,0],centered_data[:,1], yerr=centered_data[:,2]*error_scale, **H_kwargs)

    if date_unit=='days':
        plt.errorbar(centered_data[:,0],centered_data[:,1], yerr=centered_data[:,2]*error_scale, **H_kwargs)


    # plot outliers
    
    if loutliers:
        ts_outliers = __window_ts(data[loutliers,:],xmin=xmin,xmax=xmax)
        if ts_outliers==[]:return()
        
        ts_outliers[:,1]=ts_outliers[:,1]-yaxis_shift
        
        plt.errorbar(ts_outliers[:,0],ts_outliers[:,1], yerr=ts_outliers[:,2]*error_scale, fmt='o', markersize=H_kwargs['markersize'], color='r',ecolor='grey',capsize=0,linewidth=1.0)
    
    plt.grid(True)
    
    # now xmin & xmax are fixed
    
    xmin=centered_data[0,0]
    xmax=centered_data[-1,0]

    # fancy calendar 
    if date_unit=='cal':
        import matplotlib.dates as mdates
my_ts = ts.AM01.extract_periods([2016.1,2016.8])
f,ax=plt.subplots(1)
myFmt = mdates.DateFormatter('%m-%d')
ax.xaxis.set_major_formatter(myFmt)


    return([xmin,xmax])

###############################################################################
def plot(self,\
         title=None,\
         loffset=True,\
         loutliers=True,\
         verbose=False,\
         date=[],\
         yaxis=None, \
         yupaxis=None, \
         lcomponent=['N','E','U'],\
         error_scale=1.0, \
         lperiod=[[]],\
         lvline=[],\
         save_dir_plots='.', \
         save=None,\
         show=True,\
         unit='mm',\
         date_unit='decyear', \
         date_ref=0.0, center=True,\
         superimposed=None, \
         superimposed_equal_date=None, \
         plot_size=None, \
         info = [], \
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
    :param superimposed_equal_date: date at which the master and superimposed gts will be equal (default=None). date can also be a list with [date,offset_north,offset_east,offset_up]
    :param plot_size: plot size as a tuple. Default, best guess.
    :param info: info to appear in time series subtitles. A list of keywords among: 'wrms','bias','vel','period'
    :param **kwargs: any argument to be passed to plot_date (unit_date='cal') or errorbar (unit_date='decyear' or 'days')

    :note: The list of kwargs are:
      H_kwargs_plot_date={\\
                         'marker' : marker_main_symbol,\\
                         'markersize' : marker_main_size, \\
                         'markerfacecolor' : marker_main_color,\\
                         'markeredgecolor' : marker_main_color,\\
                         'markeredgewidth' : marker_main_colorlw \\
                         }
     
    H_kwargs_plot_date_superimposed={ \\
                        'color' : line_2_color,\\
                        'linewidth' : line_2_width,\\
                        'linestyle' : line_2_style,\\
                        'marker' : marker_2_symbol,\\
                        'markersize' : marker_2_size, \\
                        'markerfacecolor' : marker_2_color,\\
                        'markeredgecolor' : marker_2_color,\\
                        'markeredgewidth' : marker_2_colorlw \\
                        }
    
    
    H_kwargs_errorbar={'linewidth' : 0,\
                       'marker' : marker_main_symbol ,\
                       'markersize' : marker_main_size, \
                       'markerfacecolor' : marker_main_color,\
                       'markeredgecolor' : marker_main_color,\
                       'markeredgewidth' : marker_main_colorlw,\
                       'ecolor' : error_bar_color,\
                       'elinewidth' : error_bar_linewidth,\
                       'capsize' : error_bar_capsize }
    
    H_kwargs_errorbar_superimposed={ \
                        'color' : line_2_color,\
                        'linewidth' : line_2_width,\
                        'linestyle' : line_2_style,\
                        'marker' : marker_2_symbol,\
                        'markersize' : marker_2_size, \
                        'markerfacecolor' : marker_2_color,\
                        'markeredgecolor' : marker_2_color,\
                        'markeredgewidth' : marker_2_colorlw, \
                        'ecolor' : error_bar_color,\
                        'elinewidth' : error_bar_linewidth,\
                        'capsize' : error_bar_capsize }


    """

    # check data is not None
    from pyacs.gts.lib.errors import GtsInputDataNone
    
    try:
        if self.data is None:
            # raise exception
            raise GtsInputDataNone(inspect.stack()[0][3],__name__,self)
    except GtsInputDataNone as error:
        # print PYACS WARNING
        print( error )
        return( self )
    
    
    # cdata
    
    #if not self.cdata(data=True):
    #    print '! Cannot make plot. Problem with .data attribute.'
    #    return(self)

    # deal with **kwargs for plot settings
    
    __init_kwargs()
    
    if date_unit == 'cal':
        # **kwarg for matplotlib plot_date
        for key,value in list(kwargs.items()):
            H_kwargs_plot_date[key]=value

    if (date_unit == 'days') or (date_unit == 'decyear'):
        # **kwarg for matplotlib errorbar
        for key,value in list(kwargs.items()):
            H_kwargs_errorbar[key]=value

    # we usually want the plot in mm

    if unit=='mm':to_mm=1000.0
    if unit=='cm':to_mm=100.0
    if unit=='m':to_mm=1.0

    
    # import new modules
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    import pyacs.lib.astrotime
    
    # turn show to off by default
    plt.ioff()
    
    # some drawing default parameters
    params = {'legend.fontsize': 4}
    plt.rcParams.update(params)

    # For figure caption
    txt_component={}
    txt_component['N']='North'
    txt_component['E']='East'
    txt_component['U']='Up'

    # format for printing duration in subplot titles

    if date != []:duration=date[-1]-date[0]
    else:duration=self.data[-1,0]-self.data[0,0]
    
    # adjust tickers for the decila years case
    if date_unit=='decyear':        
        fmt=ticker.FormatStrFormatter('%6.1lf')
        if duration < 1.0:fmt=ticker.FormatStrFormatter('%6.2lf')
        if duration < 0.2:fmt=ticker.FormatStrFormatter('%6.3lf')
        if duration < 0.1:fmt=ticker.FormatStrFormatter('%6.4lf')
        
        plt.rcParams['axes.formatter.limits']=(-7,17)
    
    # ensure superimposed is either None or a list of gts
    
    if superimposed is not None:
        if not isinstance(superimposed,list):
            superimposed=[superimposed] 

        for agts in superimposed:
            if not isinstance(agts,Gts):
                raise TypeError
    
    #for subplot title information
    if self.code:stitle=self.code
    else:stitle=''

    # list of color for superimposed ts
    lcolor = ['r','g','c','m','y','k','b']


    if superimposed is not None:
        i=0
        for agts in superimposed:
            stitle+=' / '
            stitle+=agts.code+("-%s" % lcolor[i])
            i = i + 1
            
    if title is not None:
        stitle = stitle + ' ' + title

    # ensure lperiod is a list of list

    lperiod = __ensure_list_of_list(lperiod)

    # plot size

    if plot_size is None:
        if len(lcomponent)==1:plot_size=plot_size_1_component
        if len(lcomponent)==2:plot_size=plot_size_2_component
        if len(lcomponent)==3:plot_size=plot_size_3_component
    
    plt.figure(num=None, figsize=plot_size)

    # a small routine to select the data to be plotted depending on the component
    
    def __data_4_plot(gts, component, unit='mm'):
        if component=='N':data=gts.data[:,(0,1,4)]
        if component=='E':data=gts.data[:,(0,2,5)]
        if component=='U':data=gts.data[:,(0,3,6)]

        data[:,1:]=data[:,1:]*to_mm
        return(data)

    # handle superimposed_equal_date option
    
    if superimposed_equal_date is not None:

        center = False

        wgts=self.set_zero_at_date(superimposed_equal_date,offset=0.0)
        
        for i in np.arange( len(superimposed)):
            
            superimposed[i] = superimposed[i].set_zero_at_date(superimposed_equal_date,offset=0.0)
        
        wsuperimposed=superimposed
        
        
    else:
        wgts=self
        wsuperimposed=superimposed
        
    # Create frames
    
    
    for component in lcomponent:
        if component == 'N' and len(lcomponent)==3:
            nsubplot=311
        if component == 'E' and len(lcomponent)==3:
            nsubplot=312
        if component == 'U' and len(lcomponent)==3:
            nsubplot=313

        if component == 'N' and len(lcomponent)==2:
            nsubplot=211
        if component == 'E' and len(lcomponent)==2:
            nsubplot=212

        if len(lcomponent)==1:
            nsubplot=111
        
        # setup_save current subplot
        plt.subplot(nsubplot)

#        if date_unit=='decyear':        
#            plt.subplot(nsubplot).xaxis.set_major_formatter(fmt)

        if wgts.velocity is None:
            vel=None
        else:
            if component == 'N':vel=wgts.velocity[0]*to_mm
            if component == 'E':vel=wgts.velocity[1]*to_mm
            if component == 'U':vel=wgts.velocity[2]*to_mm

        str_title=__make_stitle(component,vel,__data_4_plot(wgts,component),unit, info=info)

        
        # plot the data

        if date_unit == 'cal':
                (xmin_subplot,xmax_subplot) = __plot_ts(plt,__data_4_plot(wgts,component),\
                        loutliers=wgts.outliers,\
                        date=date,\
                        date_unit=date_unit,\
                        ref_date=date_ref,\
                        yaxis=yaxis,\
                        center=center,\
                        superimposed=None,\
                        verbose=verbose,\
                        error_scale=error_scale,
                        **H_kwargs_plot_date)
    
        if (date_unit == 'days') or (date_unit == 'decyear'):
                (xmin_subplot,xmax_subplot) = __plot_ts(plt,__data_4_plot(wgts,component),\
                        loutliers=wgts.outliers,\
                        date=date,\
                        date_unit=date_unit,\
                        ref_date=date_ref,\
                        yaxis=yaxis,\
                        center=center,\
                        superimposed=None,\
                        verbose=verbose,\
                        error_scale=error_scale,
                        **H_kwargs_errorbar)

        # plot offset dates as vertical lines
    
        if loffset:

            if date_unit == 'days':
                if (date_ref is None) or (date_ref == 0):
                    # reference data is year.000
                    date_ref=float(int(wgts.data[0,0]))
                
                ldate_offsets_in_date_unit=pyacs.lib.astrotime.day_since_decyear(wgts.offsets_dates,date_ref)

            if date_unit == 'decyear':
                ldate_offsets_in_date_unit=np.array(wgts.offsets_dates)-date_ref
            
            if date_unit == 'cal':
                ldate_offsets_in_date_unit = __decyear2num(wgts.offsets_dates)


            for odate in ldate_offsets_in_date_unit:
                if  odate < xmax_subplot and  odate > xmin_subplot:
                    plt.axvline(x= odate ,linewidth=2, color='g', linestyle='--')

     
        
        # subplot title
        plt.title(str_title, x=0, horizontalalignment='left', fontsize= '11')
        
        # superimposed option

        if superimposed is not None:

            
            for i in np.arange(len(wsuperimposed)):

                agts = wsuperimposed[i]
                # color for agts

                H_kwargs_plot_date_superimposed['color'] = lcolor[i]
                H_kwargs_errorbar_superimposed['color'] = lcolor[i]
                

                # calendar date
            
                if date_unit == 'cal':
                        (xmin_subplot,xmax_subplot) = __plot_ts(plt,__data_4_plot( agts , component),\
                                loutliers = agts.outliers,\
                                date=date,\
                                date_unit=date_unit,\
                                ref_date=date_ref,\
                                yaxis=yaxis,\
                                center=center,\
                                superimposed=None,\
                                verbose=verbose,\
                                error_scale=error_scale,
                                **H_kwargs_plot_date_superimposed)

                # date_init is 'days' or 'decyear'
                
                if (date_unit == 'days') or (date_unit == 'decyear'):
                        (xmin_subplot,xmax_subplot) = __plot_ts(plt,__data_4_plot( agts ,component),\
                                loutliers = agts.outliers,\
                                date=date,\
                                date_unit=date_unit,\
                                ref_date=date_ref,\
                                yaxis=yaxis,\
                                center=center,\
                                superimposed=None,\
                                verbose=verbose,\
                                error_scale=error_scale,
                                **H_kwargs_errorbar_superimposed)

        # lperiod option: periods to be highlighted 
        
        if lperiod != [[]]:
            if date_unit == 'days':
                if (date_ref is None) or (date_ref == 0):
                    # reference data is year.000
                    date_ref=float(int(wgts.data[0,0]))

            __plot_lperiod(plt,lperiod,date_unit,date_ref)

        # vertical lines option
        if lvline != []:

            if date_unit == 'days':
                if (date_ref is None) or (date_ref == 0):
                    # reference data is year.000
                    date_ref=float(int(wgts.data[0,0]))
                
                ldate_vline_in_date_unit=pyacs.lib.astrotime.day_since_decyear(lvline,date_ref)

            if date_unit == 'decyear':
                ldate_vline_in_date_unit=np.array(lvline)-date_ref
            
            if date_unit == 'cal':
                ldate_vline_in_date_unit = __decyear2num(lvline)


            for odate in ldate_vline_in_date_unit:
                if  odate < xmax_subplot and  odate > xmin_subplot:
                    plt.axvline(x= odate ,linewidth=1, color='k', linestyle='-')

        
        # Y-label    
        plt.ylabel(unit)

        # Xticks
        if date_unit == 'decyear':
            plt.gca().xaxis.set_major_formatter(fmt)


    # X Label
    
    if date_unit=='decyear':
        if (date_ref is None) or (date_ref == 0.0):
            str_xlabel='year'
        else:
            str_xlabel=("decimal year since %s" %pyacs.lib.astrotime.decyear2datetime(date_ref))
    

    if date_unit=='days':
        
        if (date_ref is None) or (date_ref == 0):
            # reference data is year.000
            date_ref=float(int(wgts.data[0,0]))
        
        str_xlabel=("days since %s" %pyacs.lib.astrotime.decyear2datetime(date_ref))

    if date_unit=='cal':
        str_xlabel=("calendar date")

    
    plt.xlabel(str_xlabel)
    plt.suptitle(stitle)

    
    # tight_layout

    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=hspace)
    #plt.tight_layout()
    
    # end of subplots

    
    # save file as png if save option is set
    
    if save is not None:

        # import
                
        import os
        import os.path
        
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
        
        plt.savefig(fname, dpi=150, facecolor='w', edgecolor='w', orientation='portrait', papertype='a4', format='png',transparent=False, pad_inches=0.1)
    
    # show
    if show: 
        # do not block program execution when displaying figures
        plt.ion()
        plt.show()

    # returns the Gts for pyacs convention compatibility        
    new_Gts=self.copy()
    
    return(new_Gts)
    
