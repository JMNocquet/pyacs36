"""
This module gathers a few useful date functions, which can be used to read different time series formats.
pyacs.pygts uses decimal year as time defaults. Seconds for High-Rate GPS solution can also be used.
"""

import numpy as np
import datetime

def np_yyyy_mm_dd_hh_mm_ss_2_decyear(data):
    """
    converts a numpy array including year month mday hour minute sec to decimal year
    returns a 1D array
    """
    
    from pyacs.lib import  astrotime as AstroTime

    ndates=data.shape[0]

    np_dec_year=np.zeros(ndates)
    for i in range(ndates):
        uts=AstroTime.hmsmicros2uts(data[i,3], data[i,4], s=int(data[i,5]), microsecond=(data[i,5]-int(data[i,5]))/1.E6)
        np_dec_year[i]=AstroTime.cal2decyear(data[i,2], data[i,1], data[i,0], ut=0.5, uts=uts)
        
    return(np_dec_year)
    
def np_yyyy_mm_dd_hh_mm_ss_2_datetime(data):
    """
    converts a numpy array including year month mday hour minute sec to an array of python datetime.datetime object
    returns a hash
    """

    import datetime
    ndates=data.shape[0]

    np_datetime={}
    for i in range(ndates):
        # round error in second time in track e.g. -0.000001
        if data[i,5]<0 and data[i,5]>-0.0001:data[i,5]=0.0
        
        microsecond=int((data[i,5]-int(data[i,5]))*1.E6)
        np_datetime[i]=datetime.datetime(int(data[i,0]), int(data[i,1]), int(data[i,2]), int(data[i,3]), int(data[i,4]),int(data[i,5]),microsecond)
    
    return(np_datetime)

def np_datetime_2_eq_time(data,leap_sec=0.0,eq_time=0.0):
    """
    takes a hash of python datetime.datetime object and return a numpy array of seconds with respect to eq_time
    if the input array is in GPS time, providing leap_sec correct for the GPS_time - UTC delta

    :param leap_sec: number of seconds between GPS_time - UTC delta (leap_sec=17 that is GPS is ahead of UTC by 17 seconds on 13/02/2016)
    :param eq_time:  time of earthquake as a python datetime.datetime object (in UTC)
    """
    
    from datetime import timedelta
    
    corrected_data=np.zeros((len(data)))
    
    for i in sorted(data.keys()):
        #corrected_data[i]=(data[i]-timedelta(seconds=leap_sec)-eq_time).total_seconds()
        corrected_data[i]=(data[i]-eq_time).total_seconds()-leap_sec
            
    return(corrected_data)
 
def decyear2days(self,ref_date='',in_place=False):
    """
    Converts the dates of a time series from decimal years to days after a reference date
    ref_date is read by guess_date
    """

    from pyacs.lib.astrotime import guess_date
    
    if ref_date=='':
        ref_date=self.data[0,0]
    
    def np_decyear_2_days(data,ref_date):
        """
        converts a 1-D numpy array including decimal year to a 1-D numpy array of days after a reference date
        ref_date is read by guess_date
        returns a 1-D numpy array
        """
    
        from pyacs.lib import  astrotime as AT
    
        lmjd=np.array(list(map(AT.decyear2mjd,data)))
        ref_mjd=AT.decyear2mjd(guess_date(ref_date))
        
        return(lmjd-ref_mjd)

    new_Gts=self.copy()
    
    new_Gts.data[:,0]=np_decyear_2_days(self.data[:,0],ref_date)
    
    if in_place:
        self.data=new_Gts.data
    return(new_Gts)


