###################################################################
def to_displacement(self,verbose=True,base_name='vel',wdir='.',up=False):
###################################################################
    """
    print displacements every dates as gmt psvelo files 
    """
     
    # import
     
    import pyacs.lib.astrotime as AT
    import pyacs.vel_field as Velocity_Field
    import pyacs.lib.gmtpoint as GMT_Point
    import numpy as np

    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG


    import inspect

    VERBOSE("Running Sgts.%s" % inspect.currentframe().f_code.co_name)

    # get the list of dates
     
    ldate=[]
    for gts in self.lGts():
        if gts.data is None:
            WARNING("No data for time series %s" % gts.code)
            continue
        ldate=ldate+gts.data[:,0].tolist()
    
    np_dates=np.array(sorted(list(set(ldate))))
     
    VERBOSE("Found %d dates" % np_dates.shape[0])
    
    # creates the gmt files
     
    date_ref=np_dates[0]
     
    for i in np.arange(np_dates.shape[0]):
        date=np_dates[i]
    #        for date in sorted(ldate):
        file_name=wdir+'/'+("%04d_" % i)+base_name+'.gmt'
        vel=Velocity_Field.Velocity_Field(file_name=file_name)
        for gts in self.lGts():
            M=GMT_Point.GMT_Point(code=gts.code)
             
            tmp_gts=gts.extract_dates([date])
             
            if tmp_gts.data is not None:
                if tmp_gts.data.shape[0]==1:
                    M.lon=gts.lon
                    M.lat=gts.lat
                    if up:
                        M.Ve=0.0
                        M.Vn=tmp_gts.data[0,3]*1.E3
                        M.SVe=0.0
                        M.SVn=tmp_gts.data[0,6]*1.E3
                        M.SVen=0.0
                    else:
                        M.Ve=tmp_gts.data[0,2]*1.E3
                        M.Vn=tmp_gts.data[0,1]*1.E3
                        M.SVe=tmp_gts.data[0,5]*1.E3
                        M.SVn=tmp_gts.data[0,4]*1.E3
                        M.SVen=0.0
                     
                    vel.add_point(M)
         
        (mday,month,ut)=AT.decyear2cal(date)
        (noday,ut)=AT.decyear2dayno(date)    
        date_info=("step #%04d date %10.5lf %02d-%02d-%04d-%.1lf doy %03d | days since reference date +%4.1lf " % \
                   (i,date,mday,month,int(date),ut,noday,(date-date_ref)*365.25))
    
        vel.write(comment=date_info)
