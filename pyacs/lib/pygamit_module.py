"""
    module for Gamit
"""


def read_qfile(qfile):
    """Reads a qfile and returns the following information:
       - number of sites included in the analysis
       - list of sites that do not have solution
       - list of sites with large adjustment
       - total number of ambiguities
       - total number of fixed ambiguities
       - the pre/postfit values
       - Normal end in solve
    """
    
    # variables intitialization
    
    lsite_not_calculated=[]
    normal_stop=False
    lfit=[]
    
    pattern_total_ambiguity='Phase ambiguities in solution'
    pattern_NL_ambiguity_fixed='NL ambiguities resolved'
    pattern_WL_ambiguity_fixed='WL ambiguities resolved by AUTCLN'
    # TO BE FIXED
    #pattern_not_calculated='GEOC LAT  dms'
    pattern_not_calculated='POLOPOPLOPOL'
    pattern_fit=' Prefit nrms:'
    pattern_normal_stop='Normal stop in SOLVE'
    
    # open file
    
    fqfile=open(qfile)
    
    for line in fqfile:
        
        if len(line)<5:continue
        if pattern_total_ambiguity in line:n_total_ambiguity=int(line.split()[0])
        if pattern_NL_ambiguity_fixed in line:n_NL_fixed_ambiguity=int(line.split()[0])
        if pattern_WL_ambiguity_fixed in line:n_WL_fixed_ambiguity=int(line.split()[0])
        if pattern_not_calculated in line:
            if line[5]!='*':
                if line[6:10] not in lsite_not_calculated:lsite_not_calculated.append(line[6:10])
        if pattern_normal_stop:normal_stop=True
        if pattern_fit in line:
            lfit.append((float(line.split()[2]),float(line.split()[5])))
    fqfile.close()
    return(n_total_ambiguity,n_NL_fixed_ambiguity,n_WL_fixed_ambiguity,lsite_not_calculated,lfit,normal_stop)
        

def get_nearest_from_apr(aprfile,n,linclude,lexclude,X,Y,Z,dec_year):
    """Returns the n nearest site from X,Y,Z present in linclude and apr
    """
    import math
    H_distance={}
    H_site={}
    lreturn=[]
    
    
    fapr=open(aprfile)
    
    for line in fapr:
        
        if line[0]!=' ':continue
        site=line[1:5].lower()
        if site not in linclude:continue
        if site in lexclude:continue
        lline=line.split()
        (X1,Y1,Z1,_dec_year_current)=(float(lline[1]),float(lline[2]),float(lline[3]),float(lline[7]))
        distance=math.sqrt((X-X1)**2+(Y-Y1)**2+(Z-Z1)**2)
        
        if site in H_site:continue
#            dec_year1=float(H_site[site].split()[7])
#            if math.fabs(dec_year1-dec_year)<math.fabs(dec_year_current-dec_year):
#                continue
        H_site[site]=line.rstrip()
        H_distance[distance]=line.rstrip()
    
    fapr.close()
    i=0
    for distance in sorted(H_distance.keys()):
        lreturn.append(H_distance[distance])
        i=i+1
        if i>=n:break
    
    return(lreturn) 
            

def substitute_site(aprfile,site,X,Y,Z,dec_year):
    """Substitute an entry in apr file
    """
    Lline=[]
    new_line_in_apr=(" %s_GPS %14.4lf %14.4lf %14.4lf  0.00000    0.00000    0.00000 %8.4lf\n" % (site.upper(),X,Y,Z,dec_year))
    fapr=open(aprfile)
    
    for line in fapr:
        
        if line[0]!=' ':
            Lline.append(line)
            continue
        site_apr=line[1:5].lower()
        if site_apr != site.lower():
            Lline.append(line)
    fapr.close()
    Lline.append("# This line has been added automatically by pyacs\n")
    Lline.append(new_line_in_apr)
    
    fapr=open(aprfile,'w')
    fapr.writelines(Lline)
    fapr.close()
    
    
    
def update_apr(input_apr,target_apr, threshold, action):
    """ Update an apr file using another apr
        threshold is a minimum value for updating
        if action is False only check
    """
        

    # Reading the two apr
    
    f=open(input_apr)
    l_input_apr=f.readlines()
    f.close()

    f=open(target_apr)
    l_target_apr=f.readlines()
    f.close()

    # create a hash for target_apr
    H_target_apr={}
    for line in l_target_apr:
        if line[0]!=' ':continue
        site=line[1:5].upper()
        if site in H_target_apr:H_target_apr[site].append(line)
        else:H_target_apr[site]=[line]
    
    # comment to be added at the end of a line
    
    import time
    import math
    (year,month,day,_hour,_min,_sec,_weekday,_yearday,_dlsflag)=time.localtime()
    message="Added_by_pyacs_"+("%04d_%02d_%02d" % (year,month,day))
    
    # loop on input input apr lines
    lnew_apr=[]
    H_update={}
    
    for line_input in l_input_apr:
        if line_input[0]!=' ':continue
        
        
        lline_input=line_input.split()
        site_input=line_input[1:5].upper()
        (X1,Y1,Z1,dec_year_1)=(float(lline_input[1]),float(lline_input[2]),float(lline_input[3]),float(lline_input[7]))
        OK=True
        for line in H_target_apr[site_input]:
            site_apr=line[1:5].upper()
            if site_apr!=site_input:continue
            lline=line.split()
            VX2=0
            VY2=0
            VZ2=0
            if len(lline)>7:
                (X2,Y2,Z2,VX2,VY2,VZ2,dec_year_2)=\
                (float(lline[1]),float(lline[2]),float(lline[3]),float(lline[4]),float(lline[5]),float(lline[6]),float(lline[7]))
            else:
                (X2,Y2,Z2,dec_year_2)=(float(lline[1]),float(lline[2]),float(lline[3]),float(lline[7]))

            # check whether update is needed
            X2t=X2+VX2*(dec_year_1-dec_year_2)
            Y2t=Y2+VY2*(dec_year_1-dec_year_2)
            Z2t=Z2+VZ2*(dec_year_1-dec_year_2)
            distance=math.sqrt((X2t-X1)**2+(Y2t-Y1)**2+(Z2t-Z1)**2)
            
            OK=True
            if distance < threshold:OK=False

        if OK:                
            print(("=> site to be updated %s distance :%10.2lf (m)"  % (site_input.upper(), distance)))
            H_update[site_input.upper()]=line_input.rstrip()+' '+message+"\n" 

    # now doing the update

    
    for line in sorted(l_target_apr):
        site=line[1:5].upper()
        if site in H_update:
            lnew_apr.append(H_update[site])
            date_target_apr=float(line.split()[7])
            date_input_apr=float(H_update[site].split()[7])
            if math.fabs(date_target_apr-date_input_apr)>(1./365.25):lnew_apr.append(line)
            del(H_update[site])
        else:
            lnew_apr.append(line)

    # reformat entries
    
    lformatted=[]
    for line in sorted(lnew_apr):
        line.rstrip()
        lline=line.split()
        if len(lline)>=8 and line[0]==' ':
            (X,Y,Z,VX,VY,VZ,date)=list(map(float,lline[1:8]))
            fmt_line=(" %s %14.4lf %14.4lf %14.4lf  %9.5lf %9.5lf %9.5lf %8.3lf " % (lline[0],X,Y,Z,VX,VY,VZ,date))
            if len(lline)>8:
                for i in range(9,len(lline)):
                    fmt_line=fmt_line+' '+lline[i]
            
            fmt_line=fmt_line+"\n"
            lformatted.append(fmt_line)

        else:
            temp=line
            if len(temp.strip())>0:
                lformatted.append(line+"\n")
    if action=='update':
        f=open(target_apr,'w')
        f.writelines(lformatted)
        f.close()


def make_apr_doy(input_apr,output_apr, lsite, doy, year):
    """ Creates an apr file for a given list of sites and select the best coordinates according to dates
    """

    import math
    
    dec_year=year+doy/365.25

    LSITE=[]
    for site in sorted(lsite):
        LSITE.append(site.upper())
    lsite=LSITE
    
    # Reading the two apr
    
    f=open(input_apr)
    l_input_apr=f.readlines()
    f.close()

    # go through the apr
    H_kept={}
    for line_input in l_input_apr:
        if line_input[0]!=' ':continue
        lline_input=line_input.split()
        
        site=line_input[1:5].upper()
        if site not in lsite:continue
        dec_year1=float(lline_input[7])
        if site in H_kept:
            fields=H_kept[site].split()
            dec_year2=float(fields[7])
            if (math.fabs(dec_year-dec_year1)<math.fabs(dec_year-dec_year2)):
                H_kept[site]=line_input
        else:
            H_kept[site]=line_input
    
    # save results
    
    f=open(output_apr,'w')
    f.writelines(sorted(H_kept.values()))
    f.close()
        
    