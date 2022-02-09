"""
Reads GAMIT/GLOBK TRACK kinematics time series
"""
###################################################################
def read_track_NEU(self,tsdir='.',tsfile=None, leap_sec=0.0, eq_time=None, verbose=False):
###################################################################
    """
    Read a GAMIT/GLOBK Track output file generated with the option out_type NEU
    in this case dates are seconds
    by default the seconds are with respect to the first epoch of measurements
    If option leap_sec is provided with a value > 0.0, then GPS time is corrected for the difference between GPTS time and UTC 
    If eq_time is provided, it is assumed to be UTC. Expected format is YYYY:MM:MD:HH:MM:SS.S
    """

    import numpy as np
    import os

    
    # try to guess name
    if isinstance(tsfile,str):
        track_file=tsfile
    else:
        from glob import glob
        track_file=glob(tsdir+'/'+'*.NEU.'+self.code.lower()+'.LC')[0]

    if not os.path.isfile(track_file):
        print("!!! Could not open ",track_file)
        return()
    
    try:
        data=np.genfromtxt(track_file, comments='*', skip_header=2, usecols=list(range(12)))
    except:
        print("!!! Error while reading ",track_file)
        return()
    
    # handling dates
    from pyacs.gts.lib import gts_dates
    
    data_dates=data[:,:6]
    data_in_datetime=gts_dates.np_yyyy_mm_dd_hh_mm_ss_2_datetime(data_dates)
    
    #print data_in_datetime
    import datetime
    if isinstance(eq_time,str):
        [year, month, day, hour, minute, second]=list(map(int,eq_time.split(':')))
        eq_time_datetime=datetime.datetime(year, month, day, hour, minute, second)
    else:
        eq_time_datetime=data_in_datetime[0]
        
    data_in_eq_time=gts_dates.np_datetime_2_eq_time(data_in_datetime, leap_sec=leap_sec, eq_time=eq_time_datetime)
    
    data_kts=np.zeros((data.shape[0],7))
    data_kts[:,0]=data_in_eq_time

    data_kts[:,1]=data[:,6]
    data_kts[:,4]=data[:,7]

    data_kts[:,2]=data[:,8]
    data_kts[:,5]=data[:,9]
    
    data_kts[:,3]=data[:,10]
    data_kts[:,6]=data[:,11]
    
    self.data=data_kts    
