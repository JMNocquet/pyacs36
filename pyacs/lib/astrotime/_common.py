"""
Common constants and private utility functions for the astrotime package.
"""

days = [31,28,31,30,31,30,31,31,30,31,30,31]

###############################################################################
# PRIVATE ARGUMENT CHECK FUNCTIONS
###############################################################################

# Is the dayno valid?
def __daynoOK(dayno,year):
    from .leap_year import leap_year
    
    if (dayno<1 or dayno>366 or (dayno>365 and not leap_year(year))):
        raise ValueError( ('!!! doy out of range: (doy/year)=(%d,%d) ' % ( dayno,year ) ) )

# Is the month valid?
def __monthOK(month):
    import numpy as np 
    if (not isinstance(month,int)) and (not isinstance(month,np.int64)):
        raise TypeError( ('!!! bad type for month. Must be integer month= ' , month ) )
        

    if (month > 12 or month < 1):
        raise ValueError( ('!!! month out of range. Must be in the range 1-12 month=%d ' % month ) )
        

# IS the day of month OK? (assumes month IS ok - should be checked first)
def __dayOK (day,month,year):
    from .leap_year import leap_year
    
    month=month-1             # For array indexing
    if (leap_year(year)):days[1] = 29
    else:days[1] = 28

    if (day < 1 or day > days[month]):
        raise ValueError( ('!!! day of month out of range: (day/month/year)=(%d,%d,%d) ' % ( day,month,year ) ) )

# Is the day fraction OK?
def __utOK (ut):
    if (ut < 0.0 or ut >= 1.0):
        raise ValueError( ('!!! fraction day ut out of range [0;1[: ut=%.9lf  ' % ( ut ) ) )

# Is the number of seconds  OK?
def __utsOK (uts):
    if (uts < 0.0 or uts >= 24.0*60.0*60):
        raise ValueError( ('!!! numer of seconds out of range [0;24.0*60.0*60[: uts=%.9lf  ' % ( uts ) ) )
