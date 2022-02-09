from datetime import datetime,timedelta

class TimePeriodUndefinedError(Exception): pass
    

"""
Time period representation and manipulation.
From F. Peix - Geoazur (15/11/2010)
Functions added by J.-M. Nocquet - Geoazur (May 2012)
"""
class TimePeriod:

    """The beginning of this period as datetime object"""
    __start_time = None
    """The end of this period as datetime object"""
    __end_time = None 
    

    def __init__ (self,start_time=None,end_time=None):
        self.__start_time = start_time
        self.__end_time = end_time

    """
    Say if this time period is defined.

    return : True if this time period is defined and false otherwise.
    """

    def isdefined(self):
        return self.__start_time != None and self.__end_time != None

    """
    Give the start time of this time period.

    return : The beginning of this time period.
    """
    def begin(self):
        if not self.isdefined():
            raise TimePeriodUndefinedError("This time period is undefined")
        return self.__start_time

    """
    Give the start time of this time period as timestamp (epoch time).

    return : The beginning of this time period.
    """
    def epoch_begin(self):
        if not self.isdefined():
            raise TimePeriodUndefinedError("This time period is undefined")
        return int(self.begin().strftime("%s"))


    """
    Give the end time of this time period.

    return : The ending of this time period.
    """
    def end(self):
        if not self.isdefined():
            raise TimePeriodUndefinedError("This time period is undefined")
        return self.__end_time

    """
    Give the end time of this time period as timestamp (epoch time).

    return : The ending of this time period.
    """
    def epoch_end(self):
        if not self.isdefined():
            raise TimePeriodUndefinedError("This time period is undefined")
        return int(self.end().strftime("%s"))
        
    """
    Compute time intersection between this time period and a given one.

    period : The given time period 


    return : The time period representing intersection 
             between this time period and the given one.
    """
    def intersection(self,period):
        # If one period is not defined
        if not self.isdefined() or not period.isdefined():
            return TimePeriod()

        # If period don't intersect
        if  self.__end_time < period.begin() or period.end() < self.__start_time:
            return TimePeriod()
    
        # For simplicity self must begin first 
        if self.__start_time > period.begin():
            # self must begin first
            return period.intersection(self)
    
        if self.__end_time < period.end():
            return TimePeriod(period.begin(),self.__end_time)
        else:
            return TimePeriod(period.begin(),period.end())

    def has_in(self,date):
    
        """
        Tells whether a given date is in the period
        
        return: True if the given date in within the period, False otherwise
        """
        if not self.isdefined():
            return TimePeriod()
    
        if self.__start_time <= date and self.__end_time >= date:
            return True
        else:
            return False

    def display(self):
        str1=datetime.isoformat(self.__start_time)
        str2=datetime.isoformat(self.__end_time)
        return("%s -> %s" % (str1,str2)) 
        
    def get_info(self):
        str1=datetime.isoformat(self.__start_time)
        str2=datetime.isoformat(self.__end_time)
        str= ("%s -> %s" % (str1,str2)) 
        return(str)



