###################################################################
def read_gts_conf(self,gts_conf_file,verbose=False):
###################################################################
    """
    Reads a gts_conf_file
    implemented commands in the file are:
    #todo add_break      [site]  [date] # date is either [decyear] [doy year] [mday month year]
    apply_offset   [site]  [offset_north,offset_east,offset_up] [date] # offset  applied is in mm, date is either [decyear] [doy year] [mday month year]
    remove_offset   [site]  [offset_north,offset_east,offset_up] [date] # offset removed is in mm, date is either [decyear] [doy year] [mday month year]
    #todo extract_periods [site] [date1,date2] # date is either [decyear] [doy year] [mday month year]
    #todo exclude_periods [site] [date1,date2] # date is either [decyear] [doy year] [mday month year]
    #todo remove_day    [site] [date]
    """
    # import
    import numpy as np
    from pyacs.lib.astrotime import guess_date
    import pyacs.lib.astrotime as at

    New_Sgts=self.copy()
    
    conf=open(gts_conf_file,'r')
    
    H_apply_offsets={}
    H_offsets_dates={}
     
    for line in conf:
        if len(line)<5 or line[0]=='#':continue
        lline=line.split()
        # apply_offset
        if lline[0] in ['apply_offset','remove_offset']:
            (code,sdn,sde,sdu,sdate)=lline[1:]
            if lline[0] == 'apply_offset':
                dn=float(sdn)/1000.0
                de=float(sde)/1000.0
                du=float(sdu)/1000.0
    
            if lline[0] == 'remove_offset':
                dn = -float(sdn)/1000.0
                de = -float(sde)/1000.0
                du = -float(sdu)/1000.0
    
             
            date=guess_date(sdate)

            if code in H_apply_offsets:
                if verbose:
                    str_date = at.decyear2datetime(date).isoformat().split('.')[0]
                    print("-- read offset to %s at %s: N %10.2lf E %10.2lf U %10.2lf " % (code,str_date,dn*1E3,de*1E3,du*1E3))

                H_apply_offsets[code]=np.vstack((H_apply_offsets[code],np.array([date,dn,de,du])))
                H_offsets_dates[code].append(date)
            else:
                if verbose:
                    str_date = at.decyear2datetime(date).isoformat().split('.')[0]
                    print("-- read offset to %s at %s: N %10.2lf E %10.2lf U %10.2lf " % (code,str_date,dn*1E3,de*1E3,du*1E3))

                H_apply_offsets[code]=np.array([[date,dn,de,du]])
                H_offsets_dates[code]=[date]
     
    # apply commands
    for gts in self.lGts():
        if verbose:print("-- Processing ",gts.code)
        if gts.code in H_apply_offsets:
            #gts.offsets_dates+=H_offsets_dates[gts.code]
            try:
                New_Sgts.__dict__[gts.code]=gts.apply_offsets(H_apply_offsets[gts.code],verbose=verbose)
            except (RuntimeError, TypeError, NameError):
                print("!!! Error applying offset to ",gts.code)
                continue
    #            else:
    #                new_gts=Gts.copy(gts)
    #            if isinstance(new_gts,Gts):
    #                New_Sgts.append(new_gts)
    #            else:
    #                print "Error processing ",gts.code, "!!! No time series created."
    
    return(New_Sgts)
