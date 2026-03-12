"""
This module gathers a few useful date functions, which can be used to read different time series formats.
pyacs.pygts uses decimal year as time defaults. Seconds for High-Rate GPS solution can also be used.
"""

import numpy as np
import datetime

def np_yyyy_mm_dd_hh_mm_ss_2_decyear(data):
    """
    Convert numpy array [year, month, mday, hour, minute, sec] to decimal year.

    Parameters
    ----------
    data : ndarray
        Array of shape (n, 6) with year, month, day, hour, minute, sec.

    Returns
    -------
    ndarray
        1D array of decimal years.
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
    Convert numpy array [year, month, mday, hour, minute, sec] to datetime objects.

    Parameters
    ----------
    data : ndarray
        Array of shape (n, 6).

    Returns
    -------
    dict
        Mapping index -> datetime.datetime.
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
    Convert dict of datetime to seconds relative to eq_time.

    If input is in GPS time, leap_sec corrects for GPS-UTC (e.g. 17 s on 2016-02-13).

    Parameters
    ----------
    data : dict
        Mapping index -> datetime.datetime.
    leap_sec : float, optional
        GPS time minus UTC in seconds (e.g. 17).
    eq_time : datetime, optional
        Earthquake time in UTC (reference for seconds).

    Returns
    -------
    ndarray
        Seconds relative to eq_time (and leap_sec correction if applied).
    """
    
    from datetime import timedelta
    
    corrected_data=np.zeros((len(data)))
    
    for i in sorted(data.keys()):
        #corrected_data[i]=(data[i]-timedelta(seconds=leap_sec)-eq_time).total_seconds()
        corrected_data[i]=(data[i]-eq_time).total_seconds()-leap_sec
            
    return(corrected_data)
 
def decyear2days(self,ref_date='',in_place=False):
    """
    Convert time series dates from decimal year to days after a reference date.

    Parameters
    ----------
    ref_date : float or str, optional
        Reference date (decimal year or parseable by guess_date). Default first data date.
    in_place : bool, optional
        If True, modify in place; otherwise return a new Gts.

    Returns
    -------
    Gts
        Time series with dates in days (or self if in_place).
    """

    from pyacs.lib.astrotime import guess_date
    
    if ref_date=='':
        ref_date=self.data[0,0]
    
    def np_decyear_2_days(data,ref_date):
        """
        Convert 1D array of decimal year to days after a reference date.

        Parameters
        ----------
        data : ndarray
            1D array of decimal years.
        ref_date : float
            Reference date in decimal year.

        Returns
        -------
        ndarray
            1D array of days after ref_date.
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


