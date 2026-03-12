
import pyacs.lib.astrotime as AstroTime
import pyacs.sol.discontinuity as Discontinuity
import datetime
import numpy as np
import subprocess
import os


###############################################################################
def glx2snx(glx, snx, dir_snx='.'):
###############################################################################
    """
    converts a GLOBK glx file into a sinex file 
    """

    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG

    # check dir_snx exists
    if not os.path.isdir(dir_snx):
        os.makedirs(dir_snx)
    
    # check output path provided in snx
    
    cmd = 'glbtosnx ' + dir_snx + ' \'\' ' + glx + ' ' + snx
    VERBOSE("Running %s" % cmd)
    subprocess.getstatusoutput(cmd)

    
###############################################################################
def snx2ssc(snx, ssc, dir_ssc):
###############################################################################
    """
    simplifies a sinex file into a ssc file 
    """

    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG


    VERBOSE("Running snx2ssc %s %s/%s" % (snx, dir_ssc, ssc) )
    # check dir_snx exists
    if not os.path.isdir(dir_ssc):
        os.makedirs(dir_ssc)

    f_in_sinex = open(snx, 'r', encoding="latin-1")
    f_out_sinex = open(dir_ssc + '/' + ssc, 'w+')
    
    flag = False
    
    for line in f_in_sinex:
        lline = line.split()
        if (lline[0] == '+SOLUTION/ESTIMATE'):flag = True
        if (lline[0] == '+SITE/ID'):flag = True
        if flag:f_out_sinex.write("%s" % line)
        if (lline[0] == '-SOLUTION/ESTIMATE'):flag = False
        if (lline[0] == '-SITE/ID'):flag = False


###############################################################################
def glred_1_day(glx, dir_prt, dir_log, apr, eq_rename=None):
###############################################################################
    """
    run glred over a day
    """

    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG

    # glx basename
    
    glx_basename = glx.split('/')[-1].split('.')[-2]
    
    # CREATES GLOBK COMMAND FILE
    
    glred_cmd = 'glred_' + glx_basename + '.cmd'
    glorg_cmd = 'glorg_' + glx_basename + '.cmd'
    
    # Creates glred
    VERBOSE("Creating %s" % glred_cmd)
    
    fglred = open(glred_cmd, 'w')
    
    if eq_rename != None:
        fglred.write("  eq_file %s\n" % eq_rename)
        
    fglred.write('  com_file @.com\n')
    fglred.write('  sol_file @.sol\n')
    fglred.write("  apr_file %s\n" % apr)
    fglred.write('  apr_neu all 1 1 1 0 0 0\n')
    fglred.write('  apr_scale 1. 0.\n')
    
    fglred.write('  apr_wob 10 10 0 0\n')
    fglred.write('  apr_ut1 10  0\n')
    
    fglred.write('  org_opt psum pbo gdl rnrp\n')
    fglred.write("  org_cmd %s\n" % glorg_cmd)
    fglred.write('  crt_opt nopr\n')
    fglred.write('  prt_opt nopr\n')
    fglred.write('  del_scra yes\n')
    
    fglred.close()
    
    # creates sites to be used for stabilization
    
    stab_file = 'stab_sites.dat'
    
    lstab_site = []
    
    VERBOSE("Creating %s" % stab_file)
    
    fstab = open(stab_file, 'w')
    fapr = open(apr, 'r', encoding="latin-1")
    for line in fapr:
        if line[0] != ' ' or len(line) < 5:continue
        site = line.split()[0][0:4].upper()
        if site not in lstab_site:
            lstab_site.append(site)
            fstab.write(" stab_site %s_@\n" % site)
    fstab.close()
    fapr.close()
    
    # glorg file
    
    VERBOSE("Creating %s " %  glorg_cmd)
    
    fglorg = open(glorg_cmd, 'w')
    
    fglorg.write(' pos_org xtran ytran ztran xrot yrot zrot scale\n')
    fglorg.write(' stab_it  5\n')
    fglorg.write(' stab_site clear\n')
    fglorg.write(' stab_min 0.005 0.003 0.005 0.003')
    fglorg.write(" source %s\n" % stab_file)
    
    fglorg.close()
    
    # now runs glred
    log = dir_log + '/log.' + glx_basename
    prt = dir_prt + '/prt.' + glx_basename
    
    lglx_file = 'lglx_' + glx_basename
    lglx = open(lglx_file, 'w')
    lglx.write(glx)
    lglx.close()
    # ,lglx,glred_cmd
    cmd = 'glred 6 ' + prt + ' ' + log + ' ' + lglx_file + ' ' + glred_cmd
    subprocess.getstatusoutput(cmd)
    os.remove(glred_cmd)
    os.remove(glorg_cmd)
    os.remove(lglx_file)
    

###############################################################################
class SSinex:
###############################################################################
    
    """
    Sinex class: reads and manipulates a Sinex solution.
    Sinex (Solution/Technique Independent Exchange format) : Sinex format handling
    Description of the Sinex format is available at:
    http://www.iers.org/documents/ac/sinex/sinex_v210_proposal.pdf 
    """
    
###############################################################################
    def __init__ (self, name=None, estimates=None, VCV=None):
###############################################################################
        
        self.name = name
        if self.name != None:
            ll = self.name.split('/')
            self.basename = ll[-1]
        
        if estimates is not None:
            self.estimates = estimates
        else:
            self.estimates = {}
            
        self.start_time = '00:000:00000'
        self.end_time = '00:000:00000'

###############################################################################
    @classmethod
    def from_sinex(cls, snx):
###############################################################################
        """
        Constructs a PYACS Sinex instance from a P. Rebishung sinex instance
        """
        
        import pyacs.lib.astrotime
        from pyacs.sol.gpoint import Gpoint
        
        name = snx.file
                
        # estimate

        estimates = {}

        for r in snx.param:
            XYZ = snx.get_xyz([r.code], pt=[r.pt], soln=None)
            M = Gpoint(X=XYZ[0, 0], Y=XYZ[0, 1], Z=XYZ[0, 2], \
                         SX=0, SY=0, SZ=0, \
                         VX=0., VY=0., VZ=0., \
                         SVX=0., SVY=0., SVZ=0., \
                         epoch=pyacs.lib.astrotime.epoch2decyear(r.tref), code=r.code, pt=r.pt, soln=int(r.soln))

            estimates[r.code, int(r.soln)] = M
  
        # VCV
        VCV = None
        
        return cls(name, estimates, VCV)
        
###############################################################################
    def clear(self):
###############################################################################
        """
        Clear the sinex
        """
        self.estimates = {}
        self.start_time = '00:000:00000'
        self.end_time = '00:000:00000'
    
###############################################################################
    def add_site(self, M):
###############################################################################
        """
        Adds a Gpoint M in the current Sinex instance
        """
        
        self.estimates[M.code, M.soln] = M

###############################################################################
    def subset_from_code(self, lcode, lsoln=None):
###############################################################################
        """
        Takes a subset from a SINEX according to a list of code and optionally a list of associated soln.
        Returns a dictionary of Gpoint with key (code,soln)
        """
        
        SUBSET = {}

        if lsoln is not None:
            for i in range(len(lcode)):
                SUBSET[lcode[i], lsoln[i]] = self.estimates[lcode[i], lsoln[i]]
                
        else: 
            for code, soln in list(self.estimates.keys()):
                if code in lcode:
                    SUBSET[code, soln] = self.estimates[code, soln]
        
        return(SUBSET)

###############################################################################
    def subset_from_keys(self, H):
###############################################################################
        """
        Takes a subset from a SINEX according to the key of a dictionary made of (code,soln)
        Returns a dictionary of Gpoint with key (code,soln)
        """
        
        SUBSET = {}

        for code, soln in list(H.keys()):
            SUBSET[code, soln] = self.estimates[code, soln]
        
        return(SUBSET)

###############################################################################
    def apply_helmert(self, T):
###############################################################################
        """
        Applies an Helmert transformation to a SINEX instance
        returns a new SINEX instance
        """


        import pyacs.message.message as MESSAGE
        import pyacs.message.verbose_message as VERBOSE
        import pyacs.message.error as ERROR
        import pyacs.message.warning as WARNING
        import pyacs.message.debug_message as DEBUG

        new_sinex = SSinex()
        
        new_sinex.name = self.name
        new_sinex.epoch = self.epoch
        
        for code, soln in list(self.estimates.keys()):
            
            VERBOSE("Applying Helmert transformation to %s %s" % (code, soln))
            
            M = self.estimates[code, soln]
            N = M.apply_helmert(T)
            new_sinex.add_site(N)

        return(new_sinex)

###############################################################################
    def write_tsxyz(self, fname, HCOV=None):
###############################################################################
        """
        Write a tsxyz file from the current Sinex instance
        """

        import pyacs.lib.glinalg

        sinex_basename = self.name.split('/')[-1]

        f = open(fname, 'w')

        for code, soln in sorted(self.estimates.keys()):
            M = self.estimates[code, soln]

            # no uncertainties            
            if HCOV is None:
                
                f.write("%s   %10.5lf   %4s %5d   %16.4lf  %16.4lf  %16.4lf %16.5lf  %16.5lf  %16.5lf %8.4f  %8.4lf  %8.4lf \n" % \
                    (sinex_basename, self.epoch, code, soln, M.X, M.Y, M.Z, 0.001, 0.001, 0.001, 0., 0., 0.,))

            # handling uncertainties

            else:
                
                COV = HCOV[code]
                
                (corr, sigma) = pyacs.lib.glinalg.cov_to_corr(COV)
                (Sx, Sy, Sz) = sigma.flatten()
                Rxy = corr[0, 1]
                Rxz = corr[0, 2]
                Ryz = corr[1, 2]

                f.write("%s   %10.5lf   %4s %5d   %16.4lf  %16.4lf  %16.4lf %16.5lf  %16.5lf  %16.5lf %8.4f  %8.4lf  %8.4lf \n" % \
                    (sinex_basename, self.epoch, code, soln, M.X, M.Y, M.Z, Sx, Sy, Sz, Rxy, Rxz, Ryz,))
                
        f.close()

###############################################################################
    def common(self, H_Gpoint, prefit=3.5, strict=True):
###############################################################################
        """
        Returns a dictionary of Gpoint common to the current Sinex object and a dictionary (code,soln) of Gpoint
        Coordinates will be the ones from the Sinex and NOT from the list of Gpoints
        commons point are points with the same code pt and soln
        """

        import pyacs.message.message as MESSAGE
        import pyacs.message.verbose_message as VERBOSE
        import pyacs.message.error as ERROR
        import pyacs.message.warning as WARNING
        import pyacs.message.debug_message as DEBUG


        H_common_Gpoint = {}
        Hdistance = {}
        
        for code, soln in list(H_Gpoint.keys()):
            if (code, soln) in list(self.estimates.keys()):
#                if code == 'CRO1':print 'OK'
                
                M = H_Gpoint[code, soln]
                N = self.estimates[code, soln]
                distance = M.xyz_distance(N)
                if (distance > prefit):
                    VERBOSE("Bad prefit (threshold value %8.3lf) for %s: %10.1lf m" % (prefit, M.code, M.xyz_distance(N)))
                    
                else:
                    Hdistance[code, soln] = distance
                    H_common_Gpoint[code, soln] = self.estimates[code, soln]
#            else:
#                if code == 'CRO1':print ' not OK'

        if strict:
            median_distance = np.median(list(Hdistance.values()))
            
            VERBOSE("Median prefit: %10.2lf mm" % (median_distance * 1.E3))
            
            for code, soln in list(Hdistance.keys()):
                if Hdistance[code, soln] > 50. * median_distance:
                    VERBOSE("Deleting %s %s bad prefit %.1lf" % (code, soln, Hdistance[code, soln] * 1.E3 ))
                    del H_common_Gpoint[code, soln]
            
        return(H_common_Gpoint)
      
###############################################################################
    def read_section_estimate(self, lexclude=[], lonly=[], discontinuity=None, rename=None, exclude=None):
###############################################################################
        """Read the ESTIMATE section of a SINEX file.

        Returns a dictionary keyed by (CODE, SOLN) with Gpoint instance values.

        Parameters
        ----------
        lexclude : list, optional
            Site codes to exclude from reading.
        lonly : list, optional
            If non-empty, only these codes are considered.
        discontinuity : object, optional
            If provided, used to populate the SOLN field (e.g. from Discontinuity/IGS soln).
        rename : dict, optional
            Keys: sinex file; values: list of ('OLD_CODE', 'NEW_CODE'). Often from Read_Conf_Pyacs.
        exclude : object, optional
            Additional exclusion filter.

        Returns
        -------
        dict
            Dictionary with (CODE, SOLN) keys and Gpoint instance values.

        Notes
        -----
        A point in SINEX is defined by CODE, PT, DOMES and SOLN; here PT and DOMES
        are ignored for versatility.
        """
        
        import pyacs.lib.astrotime
        from pyacs.sol.gpoint import Gpoint

        import pyacs.message.message as MESSAGE
        import pyacs.message.verbose_message as VERBOSE
        import pyacs.message.error as ERROR
        import pyacs.message.warning as WARNING
        import pyacs.message.debug_message as DEBUG

        # DEAL WITH RENAME IF PROVIDED
        if rename is not None:
            VERBOSE("rename info provided for sinex: %s" % self.name)
            H_rename = {}
            # Case for a CODE rename applying for all SINEX files
            if 'all' in rename:
                for (code, new_code) in rename['all']:
                    H_rename[code] = new_code
            # Case for a CODE rename applying for the current SINEX
            if self.name in list(rename.keys()):
                for (code, new_code) in rename[self.name]:
                    H_rename[code] = new_code

        # DEAL WITH EXCLUDE IF PROVIDED
        H_exclude = []
        if exclude is not None:
            VERBOSE("exclude info provided for sinex: %s" % self.name)
            # Case for a CODE rename applying for all SINEX files
            if 'all' in exclude:
                for code in exclude['all']:
                    H_exclude.append( code )
            # Case for a CODE rename applying for the current SINEX
            if self.name in list(exclude.keys()):
                for code in exclude[self.name]:
                    H_exclude.append( code )


        # FINDING SOLUTION/ESTIMATE SECTION IN SINEX
        VERBOSE("Reading ESTIMATE section of %s" % self.name)

        fs = open(self.name, 'r', encoding="latin-1")

        in_section = False
        
        for line in fs:
            # no informative line, next line please
            if (len(line.strip()) < 2 or line[0] == '*'):continue
            
            lline = line.split()

            # initialize the flag telling whether we are in the right section
            
            if (lline[0] == '+SOLUTION/ESTIMATE'):
                in_section = True
                continue
            
            if (lline[0] == '-SOLUTION/ESTIMATE'):break
            
            # if we are not in the right section, next line please
            if (not in_section):continue

            # READING SOLUTION/ESTIMATE SECTION IN SINEX

            if lline[1] == 'STAX':
                code = str(lline[2]).upper()
                pt = str(lline[3])
                soln = int(lline[4])
                epoch = str(lline[5])

                x = float(lline[8])
                sx = float(lline[9])
                
                lepoch = lline[5].split(':')
                yr = int(lepoch[0])
                doy = int(lepoch[1])
                sod = int(lepoch[2])
                year = AstroTime.yr2year(yr)

                epoch = AstroTime.epoch2decyear(lline[5])
                
            if (lline[1] == 'STAY' and str(lline[2]) == code and str(lline[3]) == pt and int(lline[4]) == soln): 
                y = float(lline[8])
                sy = float(lline[9])
            if (lline[1] == 'STAZ' and str(lline[2]) == code and str(lline[3]) == pt and int(lline[4]) == soln): 
                z = float(lline[8])
                sz = float(lline[9])

                # WE NOW HAVE ENOUGH INFORMATION TO CREATE A NEW GPOINT

                # RENAME CASE
                
                if rename:
                    if code in list(H_rename.keys()):
                        VERBOSE("Renaming %s into %s" % (code, H_rename[code]))
    
                        code = H_rename[code]

                #### IF DISCONTINUITY IS PROVIDED THEN UPDATE THE SOLN ACCORDINGLY
                # CHANGE BY JMN 02/05/2012
    
                if discontinuity is not None:
                    
                    (mday, month) = pyacs.lib.astrotime.dayno2cal(doy, year)
                    
                    hour = sod // 3600
                    left = sod % 3600
                    minute = left // 60
                    second = left % 60
    
                    mydate = datetime.datetime(year, month, mday, hour, minute, second)
    
                    solnn = discontinuity.get_soln(code, mydate)
                    if solnn != 0:soln = solnn

                # CREATES THE GPOINT INSTANCE
    
                M = Gpoint(X=x, Y=y, Z=z, \
                         SX=sx, SY=sy, SZ=sz, \
                         VX=0., VY=0., VZ=0., \
                         SVX=0., SVY=0., SVZ=0., \
                         epoch=epoch, code=code, pt=pt, soln=soln)

                self.estimates[code, soln] = M
                
            if (lline[1] == 'VELX' and str(lline[2]) == code and str(lline[3]) == pt and int(lline[4]) == soln): 
                self.estimates[code, soln].VX = float(lline[8])
                self.estimates[code, soln].SVX = float(lline[9])
            if (lline[1] == 'VELY' and str(lline[2]) == code and str(lline[3]) == pt and int(lline[4]) == soln): 
                self.estimates[code, soln].VY = float(lline[8])
                self.estimates[code, soln].SVY = float(lline[9])
            if (lline[1] == 'VELZ' and str(lline[2]) == code and str(lline[3]) == pt and int(lline[4]) == soln): 
                self.estimates[code, soln].VZ = float(lline[8])
                self.estimates[code, soln].SVZ = float(lline[9])
            
            self.epoch = epoch

        fs.close()

        # LEXCLUDE CASE
        
        for code, soln in list(self.estimates.keys()):

            # lexclude case
            if code in lexclude:
                VERBOSE("excluding %s from reading %s" % (code,self.name))
                del self.estimates[code, soln]

            # exclude case
            if code in H_exclude:
                VERBOSE("excluding %s from reading %s" % (code,self.name))
                del self.estimates[code, soln]

            # lonly case
            if ( lonly != [] ) and ( code not in lonly ):
                del self.estimates[code, soln]

#             if type=='Apr':
#                 print '-- Reading ',self.name
#     
#                 fapr = open(self.name, 'r')
#                 
#                 for line in fapr:
#                     
#                     if line[0]!=' ':continue
#                     code=line[1:5].upper()
#                     lline=line.split()
#                     (x,y,z,epoch)=(float(lline[1]),float(lline[2]),float(lline[3]),float(lline[7]))
#                     M=Gpoint(X=x,Y=y,Z=z,\
#                              SX=0.0,SY=0.0,SZ=0.0,\
#                              epoch=epoch,code=code,pt='A',soln=1)
#                     #print 'XXXX ',M.X
#                     self.estimates.append(M)
#                     M.VX=0.0
#                     M.VY=0.0
#                     M.VZ=0.0
#                 fapr.close()
#                         
#         except IOError:# as (errno, strerror):
#             print 'cannot open:', self.name
#             sys.exit()
#         
#         except ValueError:
#             print 'Value error:', self.name
#             sys.exit()

###############################################################################
    def write_sinex(self, sinex_name, Agency_Code='IRD'):

###############################################################################
        def write_header(fs):
            Mandatory_Start = '%=SNX'
            Format_Version = 1.00
            File_Agency_Code = 'PYA'
            import datetime
            now = datetime.datetime.utcnow()
            year = now.year
            month = now.month
            mday = now.day
            
            import pyacs.lib.astrotime as AstroTime
            
            Time_Sinex_Creation_Doy = AstroTime.cal2dayno (mday, month, year)
            Time_Sinex_Creation_Yr = AstroTime.year2yr(now.year)
            Time_Sinex_Creation_Seconds = now.second + now.minute * 60 + now.hour * 3600
            
            if self.start_time:Start_Time = self.start_time
            else:Start_Time = '00:000:00000'
            if self.end_time:End_Time = self.end_time
            else:End_Time = '00:000:00000'
            
            Number_Estimates = len(self.estimates)
            header = ("%5s %4.2f %3s %02d:%03d:%05d %3s %12s %12s P %05d 1 S\n" % (\
                    Mandatory_Start, Format_Version, File_Agency_Code, \
                    Time_Sinex_Creation_Yr, Time_Sinex_Creation_Doy, Time_Sinex_Creation_Seconds, \
                    Agency_Code, Start_Time, End_Time, Number_Estimates))
            fs.write(header)

        def write_section_site(fs):
            import pyacs.lib.coordinates as Coordinates
            import pyacs.lib.units as Units

            domes_default = 'XXXXXM001'
            free_des = 'No_Site_Description___'
            fs.write('+SITE/ID\n')
            for M in self.estimates:
                # print M.code,M.X,M.Y,M.Z
                (rlon, rlat, h) = Coordinates.xyz2geo(M.X, M.Y, M.Z)
                (deg_lon, mn_lon, sec_lon) = Units.radians2deg_mn_sec(rlon, range='0-360')
                (deg_lat, mn_lat, sec_lat) = Units.radians2deg_mn_sec(rlat, range='0-360')
                formatted_line = (" %4s  A %9s P %22s %03d %02d %4.1f %03d %02d %4.1f %7.1f\n" % \
                       (M.code, domes_default, free_des, deg_lon, mn_lon, sec_lon, deg_lat, mn_lat, sec_lat, h))
                fs.write(formatted_line)
            fs.write('-SITE/ID\n')

        def write_section_estimate(fs):

            def fmt(x):
                from numpy import fabs
                tmp_string = format("%21.15E" % x)
                integer = tmp_string.split('.')[0]
                exponent = int(tmp_string.split('E')[1])
                if exponent >= 0 : sign_exponent = '+'
                else:sign_exponent = '-'
                exponent = int(fabs(exponent))
                decimal = tmp_string.split('.')[1].split('E')[0][:14]
                if integer[0] == '-':
                    sign = '-'
                    ii = integer[1]
                else:
                    sign = '+'
                    ii = integer[0]
                if sign_exponent == '+':exponent = exponent + 1
                if sign_exponent == '-':exponent = exponent - 1
                
                fmt = ("%s.%s%sE%s%02d" % (sign, ii, decimal, sign_exponent, exponent))
                return(fmt)
                    
            index = 1
            fs.write('+SOLUTION/ESTIMATE\n')
            for M in self.estimates:
                epoch = AstroTime.decyear2epoch(M.epoch)
                fs.write(" %5d STAX   %4s  A %4d %12s m    2 %21s %11.5E\n" % (index, M.code, M.soln, epoch, fmt(M.X), M.SX))
                index = index + 1
                fs.write(" %5d STAY   %4s  A %4d %12s m    2 %21s %11.5E\n" % (index, M.code, M.soln, epoch, fmt(M.Y), M.SY))
                index = index + 1
                fs.write(" %5d STAZ   %4s  A %4d %12s m    2 %21s %11.5E\n" % (index, M.code, M.soln, epoch, fmt(M.Z), M.SZ))
                if M.VX:
                    index = index + 1
                    fs.write(" %5d VELX   %4s  A %4d %12s m/y  2 %21s %11.5E\n" % (index, M.code, M.soln, epoch, fmt(M.VX), M.SVX))
                    index = index + 1
                    fs.write(" %5d VELY   %4s  A %4d %12s m/y  2 %21s %11.5E\n" % (index, M.code, M.soln, epoch, fmt(M.VY), M.SVY))
                    index = index + 1
                    fs.write(" %5d VELZ   %4s  A %4d %12s m/y  2 %21s %11.5E\n" % (index, M.code, M.soln, epoch, fmt(M.VZ), M.SVZ))
                    
            fs.write('-SOLUTION/ESTIMATE\n')
        
        def write_trailer(fs):
            fs.write('%ENDSNX\n')
            
        fs = open(sinex_name, 'w')
        write_header(fs)
        write_section_site(fs)
        write_section_estimate(fs)
        write_trailer(fs)
        fs.close()

###############################################################################
# GPoint routines
###############################################################################

    def get_STA(self, code, epoch, soln=None, soln_file=None, psd_file=None, discontinuity=None):
        
        """
        Provides the predicted coordinates for a site at a given epoch. Accounts for soln and psd if provided.
        """

        import pyacs.lib.astrotime
        import pyacs.message.message as MESSAGE
        import pyacs.message.verbose_message as VERBOSE
        import pyacs.message.error as ERROR
        import pyacs.message.warning as WARNING
        import pyacs.message.debug_message as DEBUG

        # read soln file if provided
        
        if soln_file is not None:
            from pyacs.sol.discontinuity import Discontinuities
            discontinuity = Discontinuity()
            discontinuity.read_igs_discontinuity(soln_file)  
        
        # get the soln is not provided
        
        if soln is None:
            # get soln from discontinuity
            
            if (discontinuity is not None):
        
                (mday, month, ut) = pyacs.lib.astrotime.decyear2cal(epoch)
                
                sod = int(pyacs.lib.astrotime.ut2uts(ut))
                
                hour = sod // 3600
                left = sod % 3600
                minute = left // 60
                second = left % 60

                mydate = datetime.datetime(int(epoch), month, mday, hour, minute, second)

                solnn = discontinuity.get_soln(code, mydate)
                if solnn != 0:soln = solnn

                VERBOSE("soln for site %s %.5lf %s from soln file" % (code, epoch, soln))
            
            else:
                soln = 1
                VERBOSE("soln for site %s at epoch %s is %s (non info provided)" % (code, epoch, soln))
        
        # get propagated coordinates

        try:
            M = self.estimates[code, soln]
        except:
            WARNING("%s %s not present in sinex" %  (code, soln))
            return(None)

        delta_t = (epoch - M.epoch)

        VERBOSE("propagating coordinates from %.2lf to %.4lf delta_t %.2lf" % (M.epoch, epoch, delta_t))

        X_at_epoch = M.X + delta_t * M.VX
        Y_at_epoch = M.Y + delta_t * M.VY
        Z_at_epoch = M.Z + delta_t * M.VZ
        
        # adds psd contribution
        if psd_file is not None:
            from pyacs.sinex import snxutils
            import pyacs.lib.coordinates

            sepoch = pyacs.lib.astrotime.decyear2epoch(epoch)
            # Compute ENH post-seismic deformations

            from pyacs.sinex.sinex import sinex
            psd = sinex.read(psd_file)
            (denh, _senh) = snxutils.compute_psd(psd, code, sepoch)
            if np.any(denh):
                VERBOSE("psd for %s : %.1lf mm " %( code, denh * 1000.))

            # Compute XYZ post-seismic deformations
            (l, p, _he) = pyacs.lib.coordinates.xyz2geo(M.X, M.Y, M.Z)
            R = pyacs.lib.coordinates.mat_rot_local_to_general(l, p)
            dxyz = np.dot(R, denh.reshape(3, -1)).flatten()
            # Add post-seismic deformations
            X_at_epoch = X_at_epoch + dxyz[0]
            Y_at_epoch = Y_at_epoch + dxyz[1]
            Z_at_epoch = Z_at_epoch + dxyz[2]

        return([X_at_epoch, Y_at_epoch, Z_at_epoch])
      
###############################################################################
    def print_STA(self, code, soln=None):
###############################################################################
        """
        print STA values for a given code and optional soln
        """

        import pyacs.message.message as MESSAGE
        import pyacs.message.verbose_message as VERBOSE
        import pyacs.message.error as ERROR
        import pyacs.message.warning as WARNING
        import pyacs.message.debug_message as DEBUG


        if soln is not None:
            try:
                M = self.estimates[code, soln]
                VERBOSE(("%s %02d STAX %18.4lf STAY %18.4lf STAZ %18.4lf " % (M.code, M.soln, M.X, M.Y, M.Z)))
                VERBOSE(("%s %02d STAX %18.4lf STAY %18.4lf STAZ %18.4lf " % (M.code, M.soln, M.X, M.Y, M.Z)))
            except:
                WARNING("%s %s not present in %s" % (code, soln, self.name))
        else:
            try:
                for ccode, csoln in list(self.estimates.keys()):
                    if (ccode == code):
                        M = self.estimates[code, csoln]
                        VERBOSE(("%s %02d STAX %18.4lf STAY %18.4lf STAZ %18.4lf " % (M.code, M.soln, M.X, M.Y, M.Z)))
            except:
                WARNING("%s not present in %s" % (code, self.name))

###############################################################################
    def add_to_STA(self, code, soln, add_sta):
###############################################################################
        """
        add sta [dx,dx,dz] to STA values for a given code and soln
        """

        import pyacs.message.message as MESSAGE
        import pyacs.message.verbose_message as VERBOSE
        import pyacs.message.error as ERROR
        import pyacs.message.warning as WARNING
        import pyacs.message.debug_message as DEBUG


        M = self.estimates[code, soln]
        VERBOSE(("%s %02d STAX %18.4lf STAY %18.4lf STAZ %18.4lf " % (M.code, M.soln, M.X, M.Y, M.Z)))
        VERBOSE(("adding x=%10.4lf y=%10.4lf z=%10.4lf" % (add_sta[0], add_sta[1], add_sta[2])))
        M.X = M.X + add_sta[0]
        M.Y = M.X + add_sta[1]
        M.Z = M.X + add_sta[2]
        VERBOSE(("%s %02d STAX %18.4lf STAY %18.4lf STAZ %18.4lf " % (M.code, M.soln, M.X, M.Y, M.Z)))

###############################################################################
    def change_STA(self, code, soln, add_sta):
###############################################################################
        """
        change STA values for a given code and soln
        """

        import pyacs.message.message as MESSAGE
        import pyacs.message.verbose_message as VERBOSE
        import pyacs.message.error as ERROR
        import pyacs.message.warning as WARNING
        import pyacs.message.debug_message as DEBUG

        M = self.estimates[code, soln]
        VERBOSE(("%s %02d STAX %18.4lf STAY %18.4lf STAZ %18.4lf " % (M.code, M.soln, M.X, M.Y, M.Z)))
        M.X = add_sta[0]
        M.Y = add_sta[1]
        M.Z = add_sta[2]
        VERBOSE(("%s %02d STAX %18.4lf STAY %18.4lf STAZ %18.4lf " % (M.code, M.soln, M.X, M.Y, M.Z)))

###############################################################################
    def site(self, code, soln):
###############################################################################
        """
        return site for a given code and soln
        """
        return(self.estimates[code, soln])

###############################################################################
    def info_gpoint(self, code, soln=None):
###############################################################################
        """
        Print info for a given code and optionally soln
        """

        import pyacs.message.message as MESSAGE
        import pyacs.message.verbose_message as VERBOSE
        import pyacs.message.error as ERROR
        import pyacs.message.warning as WARNING
        import pyacs.message.debug_message as DEBUG
        import pyacs.debug
        from icecream import ic

        if soln is not None:
            M = self.estimates[code, soln]
            if pyacs.debug():
                ic(code, 'code, pt soln epoch', M.code, M.pt, M.soln, M.epoch)
                ic(code, 'XYZ', M.posxyz())
                ic(code, 'VXYZ', M.velxyz())
        else:
            for ccode, csoln in list(self.estimates.keys()):
                if code == ccode:
                    M = self.estimates[code, csoln]
                    if pyacs.debug():
                        ic(code, 'code, pt soln epoch', M.code, M.pt, M.soln, M.epoch)
                        ic(code, 'XYZ', M.posxyz())
                        ic(code, 'VXYZ', M.velxyz())

###############################################################################
    def lcode(self):
###############################################################################
        """
        Returns a list of all point codes in current Sinex object
        """
        lcode = []
        for M in list(self.estimates.values()):
            if (M.code not in lcode):lcode.append(M.code)
        return(lcode)

###############################################################################
    def read_apr(self, lexclude=[], discontinuity=None, rename=None):
###############################################################################
        """Read a Globk apr file.

        Parameters
        ----------
        lexclude : list, optional
            Site codes to exclude from reading.
        discontinuity : object, optional
            If provided, used to populate the SOLN field (e.g. from Discontinuity/IGS soln).
        rename : dict, optional
            Keys: sinex file; values: list of ('OLD_CODE', 'NEW_CODE'). Often from Read_Conf_Pyacs.

        Returns
        -------
        dict
            Dictionary with (CODE, SOLN) keys and Gpoint instance values.
        """
        from icecream import ic

        import pyacs.lib.astrotime
        from pyacs.sol.gpoint import Gpoint
        import pyacs.message.message as MESSAGE
        import pyacs.message.verbose_message as VERBOSE
        import pyacs.message.error as ERROR
        import pyacs.message.warning as WARNING
        import pyacs.message.debug_message as DEBUG

        # DEAL WITH RENAME IF PROVIDED
        
        if rename is not None:
            
            VERBOSE("Rename info provided for apr file: %s" % self.name)

            H_rename = {}

            # Case for a CODE rename applying for all SINEX files
            if 'all' in rename:
                
                for (code, new_code) in rename['all']:
                    H_rename[code] = new_code
            
            # Case for a CODE rename applying for the current SINEX
            
            if self.name in list(rename.keys()):

                for (code, new_code) in rename[self.name]:
                    H_rename[code] = new_code
        
        # READING APR FILE
        
        VERBOSE("Reading Globk apr file %s" % self.name)

        try:
            APR_VALUE = np.genfromtxt(self.name, comments='#', usecols=(1,2,3,4,5,6,7,8,9,10,11,12,13))
            APR_NAME  = np.genfromtxt(self.name, comments='#', usecols=(0), dtype=str)
        except:
            ERROR(("Could not read Globk format apr file: %s" % self.name), exit=True)

#        ic(APR_VALUE)
        for i in np.arange( APR_VALUE.shape[0])  :
            VERBOSE("processing %s" % APR_NAME[i][:4])
            [x,y,z,vx,vy,vz,epoch, sx,sy,sz,svx,svy,svz]= APR_VALUE[i,:]
            M=Gpoint(X=x,Y=y,Z=z,\
                     SX=sx,SY=sy,SZ=sz,\
                     VX=vx,VY=vy,VZ=vz,SVX=svx,SVY=svy,SVZ=svz, \
                     epoch=epoch,code=APR_NAME[i][:4],pt='A',soln=1)
#            ic(M.VX)
            self.estimates[ APR_NAME[i][:4], 1 ] = M
            
                 
