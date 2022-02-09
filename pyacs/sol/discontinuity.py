import datetime
import pyacs.lib.astrotime as AstroTime
import pyacs.lib.timeperiod as TimePeriod
import pyacs.sol.soln as Soln

class Discontinuities:
    
    """Reads an IGS Discontinuity file"""

    
    def __init__ (self):
        self.lsite={}
        self.lsite_offsets_dates={}
        
    def read_igs_discontinuity(self, discontinuities_file):
        
        fdiscon=open(discontinuities_file,'r')
        
        print('-- Reading ',discontinuities_file)
        ok=False
        for line in fdiscon:
            #print line
            if line[0:23] == '+SOLUTION/DISCONTINUITY':ok=True;continue
            if line[0:23] == '-SOLUTION/DISCONTINUITY':ok=False;break
            if not ok:continue
            #print line
            code=line[1:5]
            if not code in self.lsite:self.lsite[code]=Soln.Soln()
            pt=line[7:8]
            soln=int(line[9:13])
            technique=line[14:15]
            str_date_start=line[16:28]
            str_date_end=line[29:41]
            type_discon=line[42:43]
            #print (":%s:%s:%04d:%s:%s:%s:%s:\n" % (code,pt,soln,technique,str_date_start,str_date_end,type_discon))
            if type_discon=='V':continue
            # Start date
            if str_date_start == '00:000:00000' or str_date_start == '00:999:00000':
                year=1980;doy=1;sec=1;
            else:
                (yr,doy,sec)=list(map(int,str_date_start.split(':')))
                if yr>80:year=1900+yr
                else:year=2000+yr
                if sec > 86399:sec=86399
                
            (mday,month)=AstroTime.dayno2cal(doy,year)
            hour   = sec // 3600
            left   = sec % 3600
            minute = left // 60
            second = left % 60
            #print (code,year, month, mday,hour,minute,second)
            st_date=datetime.datetime(year, month, mday,hour,minute,second)
            
            if str_date_start != '00:000:00000':
                uts=24.0*hour+60.0*minute+second
                offset_date=AstroTime.dayno2decyear(doy,year,ut=AstroTime.uts2ut(uts))
                if code in self.lsite_offsets_dates:
                    if offset_date not in self.lsite_offsets_dates[code]:
                        self.lsite_offsets_dates[code].append(offset_date)
                else:
                    self.lsite_offsets_dates[code]=[offset_date]
                    
            # End date
            if str_date_end == '00:000:00000':
                year=2100;doy=1;sec=1;
            else:
                (yr,doy,sec)=list(map(int,str_date_end.split(':')))
                if yr>80:year=1900+yr
                else:year=2000+yr
                if sec > 86399:sec=86399
            (mday,month)=AstroTime.dayno2cal(doy,year)
            hour=sec // 3600
            left=sec % 3600
            minute=left // 60
            second=left % 60
            
            #print (code,year, month, mday,hour,minute,second)
            
            end_date=datetime.datetime(year, month, mday,hour,minute,second)

            if str_date_end != '00:000:00000':
                uts=24.0*hour+60.0*minute+second
                offset_date=AstroTime.dayno2decyear(doy,year,ut=AstroTime.uts2ut(uts))
                if code in self.lsite_offsets_dates:
                    if offset_date not in self.lsite_offsets_dates[code]:
                        self.lsite_offsets_dates[code].append(offset_date)
                else:
                    self.lsite_offsets_dates[code]=[offset_date]


            period=TimePeriod.TimePeriod(start_time=st_date,end_time=end_date)
            #print code,soln,period
            self.lsite[code].append_period(soln,period)
            

    def display(self, code=None):
        if not code:
            for site in sorted (self.lsite.keys()):
                #print "site: ",site
                lsoln=self.lsite[site].periods
                #print lsoln
                for soln in sorted(lsoln.keys()):
                    str=lsoln[soln].get_info()
                    print(("site:%s  soln #%s %s" % (site,soln,str))) 
        else:
            lsoln=self.lsite[code].periods
                #print lsoln
            for soln in sorted(lsoln.keys()):
                str=lsoln[soln].get_info()
                print(("site:%s  soln #%s %s" % (code,soln,str))) 
               
        
    def get_soln(self,code,date):
        """Returns the soln code corresponding to the given site and date"""
        if not code in self.lsite:return(0)
        else:
            lsoln=self.lsite[code].periods
            for soln in sorted(lsoln.keys()):
                period=lsoln[soln]
                if period.has_in(date):return(soln)
        return(0)
            
            
        
            