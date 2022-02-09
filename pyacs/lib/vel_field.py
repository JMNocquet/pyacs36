class Velocity_Field:
    
    """
    Velocity_Field class: reads a velocity field from a GMT psvelo file 
    and provides methods to manipulate velocity field
    """
    
###################################################################
    def __init__ (self, file_name=None, name=None, lgmt_points=None, verbose=False):
###################################################################

        self.name = name
        self.file_name = file_name
        if lgmt_points:
            self.sites = lgmt_points
        else:
            self.sites = []

    @classmethod        
###################################################################
    def read(cls, file_name=None, lexclude=[], lonly=[], verbose=False):
###################################################################
        
        """
        Reads a GMT psvelo file
        """

        # import
        import numpy as np
        
        # init
        
        vf = Velocity_Field()

        # fake 4-letters code generation using hexadecimal
        def __gen_fake_code__(n):
            
            FAKE = []
            for i in np.arange(n):
                fake_code = ("%4s" % hex(i).split('x')[-1].replace('L', '')).replace(' ', '0')
                FAKE.append(fake_code.upper())
            
            return(np.array(FAKE))

        # reads psvelo file

        if verbose:
            print("-- Reading GMT psvelo file: %s " % file_name)
        
        try:
            np_vel = np.array(np.mat(np.genfromtxt(file_name, comments='#')))
        except:
            raise IOError("!!! Could not read file: %s" % file_name)
        
        # empty psvelo file
        if np_vel.size == 0:
            return( vf )
        
        if (np_vel.shape[1] == 8):
            if verbose:
                print("-- file %s has 8 columns" % file_name)

            np_vel = np.delete(np_vel, -1, axis=1)
            np_code = np.array(np.mat(np.genfromtxt(file_name, comments='#', usecols=(7), dtype=str))).flatten()
            
            code =''
            np_sort = np.sort( np_code )
            for i in np.arange( np_sort.shape[0] ): 
                if np_sort[i] == code: 
                    print("!!!",code) 
                code =  np_sort[i] 

            
        elif (np_vel.shape[1] == 3):

            if verbose:
                print("-- file %s has 3 columns" % file_name)

            np_vel = np.delete(np_vel, -1, axis=1)
            np_code = np.array(np.mat(np.genfromtxt(file_name, comments='#', usecols=(2)))).flatten()

        elif (np_vel.shape[1] not in [3, 8]):
            np_code = __gen_fake_code__(np_vel.shape[0])
        else:
            raise IOError("!!! Could not decipher file content: %s", file_name)

        # populates velocity field
        
        from pyacs.lib.gmtpoint import GMT_Point

        lgmt_points = []

        for i in np.arange(np_vel.shape[0]):

            code = np_code[i]
            
            if np_vel.shape[1] >= 7:
                lon, lat, Ve, Vn, SVe, SVn, SVen = np_vel[i, :]
                M = GMT_Point(lon=lon, lat=lat, Ve=Ve, Vn=Vn, SVe=SVe, SVn=SVn, SVen=SVen, code=code)
            else:
                lon, lat = np_vel[i, :]
                M = GMT_Point(lon=lon, lat=lat, code=code)

            if verbose:
                M.get_info(display=True)
            
        # tests whether site will be added
        
            if lonly != []:
                if M.code in lonly:
                    lgmt_points.append(M)
    
            else:
                if lexclude != []:
                    if M.code not in lexclude:
                        lgmt_points.append(M)
                else:
                    lgmt_points.append(M)
    
        vf.file_name = file_name
        vf.sites = lgmt_points
         
        return vf 

###################################################################
    def info(self , details=True):
###################################################################
        """
        Print basics information
        """
        
        print("-- velocity field from: %s " % self.file_name)
        print("-- number of sites    : %d " % len(self.sites))
        
        if details:
            import numpy as np
            n_site_per_line = 10
            print('-- list of sites: ')
            A = np.append(np.array(self.lcode()), ['    '] * (n_site_per_line - np.remainder(len(self.lcode()) , n_site_per_line))).reshape(-1, n_site_per_line)
            print(np.array2string(A, separator='').replace('[', ' ').replace(']', ' ').replace('\'', ' '))
        
###################################################################

###################################################################
    def add_point(self, M):
###################################################################
        """
        Appends a GMT_Point to a velocity field
        """
        
        self.sites.append(M)

        return(self)

###################################################################
    def remove_point(self, code):
###################################################################
        """
        Removes the a GMT_Point to a velocity field given its code
        
        :param code: 4-characters code for the point to be removed
        """
        
        self = self.subset(lexclude=[code]) 
        
        return(self)

###################################################################
    def write(self, out_file=None, lexclude=[], verbose=True, comment='', up=False):
###################################################################
        """
        Writes a GMT psvelo file
        a list site name to be excluded can be provided
        """

        if (out_file is None):
            print("! No file name provided. Using tmp_vel.gmt")
            out_file = 'tmp_vel.gmt'

        if verbose:
            print("-- Writing GMT psvelo file: %s " % out_file)
        fs = open(out_file, 'a')
        
        if up:
            comment = comment + ' | UP component'

        fs.write('# ' + comment + '\n')

        if up:

            for M in self.sites:
                fs.write("%10.5lf  %10.5lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %s\n" % (M.lon, M.lat, 0.0, M.Vu, 0.0, M.SVu, M.SVen, M.code))

        else:
            
            for M in self.sites:
                fs.write("%10.5lf  %10.5lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %s\n" % (M.lon, M.lat, M.Ve, M.Vn, M.SVe, M.SVn, M.SVen, M.code))

        fs.close()
                        
###################################################################
    def nsites(self):
###################################################################
        """
        Returns the number of sites in velocity field
        """
        
        n = len(self.sites)
        return n

#    def lsite(self):
#        """Returns the site from the velocity field as a list GMT_Point object """
#        for M in self.sites:
#            if M.code == code:
#                return M

###################################################################
    def l_GMT_Point(self):
###################################################################
        
        """
        Returns the velocity field as a list of GMT_Point object 
        """
        
        lsite = []
        for M in self.sites:lsite.append(M)
        
        return lsite

###################################################################
    def print_info_site(self, code, verbose=False):
###################################################################
        """
        Print the info for a given site by his code
        """
        for M in self.sites:
            if M.code == code:
                M.get_info(legend=verbose , display=True)

###################################################################
    def subset(self, lonly=None, lexclude=None):
###################################################################

        """
        Returns a new Velocity_Field object from a subset of sites
        
        :param lonly,lexclude: list of site codes to be included or excluded
        
        
        """
        
        lGpoint = []

        for M in self.sites:
        # tests whether site will be added

            if lonly != None:
                if M.code in lonly:
                    lGpoint.append(M)
    
            else:
                if lexclude != None:
                    if M.code not in lexclude:
                        lGpoint.append(M)
        
        return Velocity_Field(lgmt_points=lGpoint)
    
###################################################################
    def radial(self, center):
###################################################################

        """
        Returns a Velocity Field whose components are radial and tangential with respect to a given center
        
        :param center: numpy array [longitude,latitude] of center
        """
        
        import numpy as np
        import copy
        
        lnew_point = []
        
        for code in self.lcode():
            
            gmt_point = self.site(code)
            
            radial_unit_vector = np.array([gmt_point.lon - center[0], gmt_point.lat - center[1]])
            radial_unit_vector = radial_unit_vector / np.sqrt(np.sum(radial_unit_vector ** 2))

            tangential_unit_vector = np.array([radial_unit_vector[1], -radial_unit_vector[0]])
            
            new_gmt_point = copy.copy(gmt_point)
            new_gmt_point.Ve = np.sum(np.array([gmt_point.Ve, gmt_point.Vn]) * radial_unit_vector)
            new_gmt_point.Vn = np.sum(np.array([gmt_point.Ve, gmt_point.Vn]) * tangential_unit_vector)
        
            lnew_point.append(new_gmt_point)

        new_vel = Velocity_Field(lgmt_points=lnew_point)
        return new_vel

###################################################################
    def lcode(self):
###################################################################
        """
        Returns a list of all point codes in current Velocity Field object
        """
        
        lcode = []
        for M in self.sites:
            if (M.code not in lcode):lcode.append(M.code)
        return(lcode)

###################################################################
    def site(self, code):
###################################################################

        """
        Returns a site in current Field object as a GMT_Point object
        """
        
        for M in self.sites:
            if (M.code == code):return(M)
        return(None)
    
###################################################################
    def calc_pole(self, plates, method='WLS' , verbose=False):
###################################################################
        """
        Performs Euler pole calculation
        
        :param plates: dictionary with the name of plate as key and list of sites for each plate as values
        :param method: choose among weighted least-squares 'WLS' or 'L1'
        
        :output: for each plate, the following files will be created:
                - euler_stat_plate_name.dat: Euler pole values and statistics
                - plate_name.gmt: velocities (GMT psvelo format) with respect to the calculated Euler pole
                - plate_name.shp: velocities (shapefile  format) with respect to the calculated Euler pole
                - eulers_sum.dat: summary of Euler pole values
        """
        
        # import 
        import pyacs.lib.shapefile
        import pyacs.lib.euler

        # loop on plates

        # dictionary of W and VCV
        H_w = {}
        
        for plate, lsite in sorted(plates.items()):
            print("-- processing plate: %s with %8d sites " % (plate, len(lsite)))
            
            if len(lsite) < 3:
                # error
                pass
            else:
                vel_plate = self.subset(lsite) 
                H_w[ plate ] = vel_plate.pole(method=method, exp=plate, log='euler_stat_' + plate + '.dat')
                
                # get the new velocity field 
                new_vel = self.substract_pole(H_w[ plate ][0], 'rot')
                
                # write it
                outgmt = plate + '.gmt'
                new_vel.write(out_file=plate + '.gmt')
                pyacs.lib.shapefile.psvelo_to_shapefile(outgmt , plate , verbose=verbose)
 
        # euler_sum
        print('-- writing euler_sum.dat')
        euler_sum_file = open('euler_sum.dat', 'w')

        euler_sum_file.write("# input file: %s\n" % self.file_name)
        euler_sum_file.write("# Euler poles\n")
        euler_sum_file.write("#-------------------------------------------------------------------------------------------------------------\n")
        
        for plate, lsite in sorted(plates.items()):

            [wx, wy, wz] = H_w[ plate ][0]
            VCV_POLE_X = H_w[ plate ][1]
            (llambda, phi, omega) = pyacs.lib.euler.rot2euler(wx, wy, wz)
            (max_sigma, min_sigma, azimuth, s_omega) = pyacs.lib.euler.euler_uncertainty([wx, wy, wz], VCV_POLE_X)

            euler_sum_file.write("%-10s | %8.3lf %8.3lf %6.3lf | %5.2lf %5.2lf %5.1lf %5.3lf | %11.5G %11.5G %11.5G\n" \
                                 % (plate, llambda, phi, omega, max_sigma, min_sigma, azimuth, s_omega , wx, wy, wz))
 
        euler_sum_file.write("#-------------------------------------------------------------------------------------------------------------\n")
        euler_sum_file.write("# Relative Euler poles\n")
        euler_sum_file.write("# Euler poles\n")
        euler_sum_file.write("#-------------------------------------------------------------------------------------------------------------\n")

        for plate_ref, _lsite_ref in sorted(plates.items()):

            [wx_ref, wy_ref, wz_ref] = H_w[ plate_ref ][0]
            VCV_POLE_X_REF = H_w[ plate_ref ][1]

            for plate, lsite in sorted(plates.items()):
    
                if plate != plate_ref:
    
                    [wx, wy, wz] = H_w[ plate ][0]
                    VCV_POLE_X = H_w[ plate ][1]

                    wx = wx - wx_ref 
                    wy = wy - wy_ref 
                    wz = wz - wz_ref 
            
                    VCV_POLE_X = VCV_POLE_X + VCV_POLE_X_REF

                    (llambda, phi, omega) = pyacs.lib.euler.rot2euler(wx, wy, wz)
                    (max_sigma, min_sigma, azimuth, s_omega) = pyacs.lib.euler.euler_uncertainty([wx, wy, wz], VCV_POLE_X)

                    euler_sum_file.write("%-5s wrt %-5s | %8.3lf %8.3lf %6.3lf | %5.2lf %5.2lf %6.1lf %5.3lf | %11.5G %11.5G %11.5G\n" \
                                         % (plate, plate_ref, llambda, phi, omega, max_sigma, min_sigma, azimuth, s_omega , wx, wy, wz))

        euler_sum_file.write("#-------------------------------------------------------------------------------------------------------------\n")
        euler_sum_file.close()
                    
###################################################################
    def pole(self, lexclude=[], method='WLS', exp='pLate', log=None):
###################################################################
        """
        Euler pole calculation; available optimization methods are WLS: weighted least squares, LS: least-squares, Dikin: L1 linear estimator
        """
        
        import numpy as np
        import pyacs.lib.coordinates as Coordinates
        import numpy.linalg as nplinalg
        import pyacs.lib.euler as euler

        np.set_printoptions(precision=2, threshold=10000, linewidth=250)
        
        # Gets dimensions for the systems
        if lexclude:
            nsites = 0
            for M in self.sites:
                if not(M.code in lexclude):
                    nsites = nsites + 1
        else:
            nsites = self.nsites()
        
        A = np.zeros([2 * nsites, 3], float)
        B = np.zeros([2 * nsites, 1], float)
        VCV_B = np.zeros([2 * nsites, 2 * nsites], float)
        
        # starts loop on sites to build the linear system

        index = 0
        H_sites = {}
        for M in self.sites:
            if not(M.code in lexclude):
                # XYZ -> local frame rotation matrix
                H_sites[M.code] = M
                R = Coordinates.mat_rot_general_to_local(np.radians(M.lon), np.radians(M.lat))
                (x, y, z) = Coordinates.geo2xyz(np.radians(M.lon), np.radians(M.lat), M.he)
                # Observation equation in local frame
                
                Ai = np.zeros([3, 3], float)
                Ai[0, 1] = z
                Ai[0, 2] = -y
                Ai[1, 0] = -z
                Ai[1, 2] = x
                Ai[2, 0] = y
                Ai[2, 1] = -x
                
                Ai = Ai / 1000.0
                
                RAi = np.dot(R, Ai)

                # fill the observation i into the general design matrix A
                
                A[2 * index, 0] = RAi[0, 0]
                A[2 * index, 1] = RAi[0, 1]
                A[2 * index, 2] = RAi[0, 2]
                A[2 * index + 1, 0] = RAi[1, 0]
                A[2 * index + 1, 1] = RAi[1, 1]
                A[2 * index + 1, 2] = RAi[1, 2]
        
                B[2 * index, 0] = M.Ve
                B[2 * index + 1, 0] = M.Vn
                
                M.assign_index(index)

                VCV_B[2 * index, 2 * index] = M.SVe ** 2
                VCV_B[2 * index + 1, 2 * index + 1] = M.SVn ** 2
                cov = M.SVe * M.SVn * M.SVen
                VCV_B[2 * index + 1, 2 * index] = cov
                VCV_B[2 * index, 2 * index + 1] = cov

                index = index + 1
        
        # solves linear system
        
        P = nplinalg.inv(VCV_B)

        ATP = np.dot((A.T), P)
        
        N = np.dot(ATP, A)
        M = np.dot(ATP, B)
        
        # Inversion
        Q = nplinalg.inv(N)
        X = np.dot(Q, M)

        # Model Prediction
        MP = np.dot(A, X)


        # Residuals
        RVen = (B - MP )

        VCV_POLE_X = Q * 1.E-12
        chi2 = np.dot(np.dot(RVen.T, P), RVen)
        dof = 2 * nsites - 3
        reduced_chi2 = np.sqrt(chi2 / float(dof))

        (wx, wy, wz) = (X[0, 0] * 1.E-6, X[1, 0] * 1.E-6, X[2, 0] * 1.E-6)

        # Convert rotation rate in the geocentric cartesian frame to Euler pole in geopgraphical coordinates
        
        (llambda, phi, omega) = euler.rot2euler(wx, wy, wz)
        (max_sigma, min_sigma, azimuth, s_omega) = euler.euler_uncertainty([wx, wy, wz], VCV_POLE_X)

        # Write log

        if (not log):flog_name = exp + '.log'
        else:flog_name = log
        flog = open(flog_name, 'w')
        
        flog.write("\nROTATION RATE VECTOR\n")
        flog.write("Wx (rad/yr): %11.5G +- %10.5G\n" % (wx, np.sqrt(VCV_POLE_X[0, 0])))
        flog.write("Wy (rad/yr): %11.5G +- %10.5G\n" % (wy, np.sqrt(VCV_POLE_X[1, 1])))
        flog.write("Wz (rad/yr): %11.5G +- %10.5G\n" % (wz, np.sqrt(VCV_POLE_X[2, 2])))
        
        # flog.write("Wx Wy Wz : %10.5G %10.5G %10.5G\n", % (wx,wy,wz))
        flog.write("ASSOCIATED VARIANCE-COVARIANCE MATRIX (rad/yr)**2\n")
        flog.write("        Wx          Wy            Wz\n")
        flog.write("---------------------------------------\n")
        flog.write("Wx |%11.5G %11.5G %11.5G\n" % (VCV_POLE_X[0, 0], VCV_POLE_X[1, 0], VCV_POLE_X[2, 0]))
        flog.write("Wy |           %11.5G %11.5G\n" % (VCV_POLE_X[1, 1], VCV_POLE_X[2, 1]))
        flog.write("Wz |                      %11.5G\n" % VCV_POLE_X[2, 2])
        flog.write("---------------------------------------\n")

        flog.write("\nEULER POLE\n")
        flog.write("longitude (dec. degree)    : %6.2lf\n" % llambda)
        flog.write("latitude  (dec. degree)    : %6.2lf\n" % phi)
        flog.write("angular velocity (deg/Myr ): %6.3lf\n" % omega)
        flog.write("ASSOCIATED ERROR ELLIPSE\n")
        flog.write("semi major axis            : %5.2lf\n" % max_sigma)
        flog.write("semi minor axis            : %5.2lf\n" % min_sigma)
        flog.write("azimuth of semi major axis : %5.1lf\n" % azimuth)
        flog.write("std(angular velocity)      : %5.3lf\n" % s_omega)

        # Calculating a few statistics

        flog.write("\nSTATISTICS\n")
        flog.write("----------\n")

        flog.write("Number of sites     = %10d\n" % nsites)
        flog.write("Chi**2              = %10.1lf\n" % chi2)
        flog.write("Reduced Chi**2      = %10.1lf\n" % reduced_chi2)
        flog.write("Deg. of. freedom    = %10d\n" % dof)
        if (dof == 0): 
            flog.write("A post. var. factor = No meaning (dof=0)\n")
        else:
            sigma_0 = np.sqrt(chi2 / dof)
            flog.write("A post. var. factor = %10.1lf\n" % sigma_0)

        # rms & wrms

        flog.write("\nRESIDUALS\n")
        flog.write("----------\n")
        flog.write("site                    pred_e     pred_n         Ve         Vn       R_ve       R_vn       S_ve       S_vn      RN_ve      RN_vn\n")
        flog.write("------------------------------------------------------------------------------------------------------------------------------------\n")

        cpt_site = 0

        rms = 0.0
        wrms = 0.0
        dwrms = 0.0
        
        for site, M in sorted(H_sites.items()):
            M = H_sites[site]
            if not(M.code in lexclude):
                i = M.get_index()

                mpe = MP[2 * i, 0]
                mpn = MP[2 * i + 1, 0]

                oe = B[2 * i, 0]
                on = B[2 * i + 1, 0]

                rve = RVen[2 * i, 0]
                rvn = RVen[2 * i + 1, 0]
                sve = M.SVe
                svn = M.SVn
 
                flog.write("%-19s %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf\n" % \
                           (M.code, mpe,mpn,oe,on,rve, rvn, sve, svn, rve / sve, rvn / svn))
            
                cpt_site = cpt_site + 1
 
                rms = rms + rve ** 2 + rvn ** 2
                wrms = wrms + (rve / sve) ** 2 + (rvn / svn) ** 2
                dwrms = dwrms + 1.0 / sve ** 2 + 1.0 / svn ** 2
        flog.write("------------------------------------------------------------------------------------------------------------------------------------\n")

        rms = np.sqrt(rms / 2 / cpt_site)
        wrms = np.sqrt(wrms / dwrms)

        flog.write("rms = %10.2lf mm/yr    wrms = %10.2lf mm/yr\n\n" % (rms, wrms))

        flog.close()

        print("-- Rotation rate results logged in %s" % flog_name)
        
        return(X * 1.E-6 , VCV_POLE_X)

###################################################################
    def substract_pole(self, W=None, type_euler='rot'):
###################################################################
        """
        substract the prediction of an Euler pole to the velocity field
        """

        import pyacs.lib.coordinates as Coordinates
        import pyacs.lib.euler as euler

        import numpy as np

        W = np.array(W).reshape(3,1)

        ltype = ['euler', 'rot']  # type is either Euler pole (lon. lat in dec. degrees and deg/Myr or cartesian rotation rate vector in rad/yr)
        if (type_euler not in ltype):
            raise ValueError(
                "!!! type_euler should be either euler or rot. type_euler=%s", type_euler)

        if (type_euler == 'euler'):
            # convert to rotation rate vector
            (lon, lat, omega) = (W[0, 0], W[1, 0], W[2, 0])
            (wx, wy, wz) = euler.euler2rot(lon, lat, omega)
            W = np.zeros([3, 1])
            W[0, 0] = wx
            W[1, 0] = wy
            W[2, 0] = wz

        # create a new velocity field

        lGpoint = []

        for M in self.sites:
                R = Coordinates.mat_rot_general_to_local(np.radians(M.lon), np.radians(M.lat))
                (x, y, z) = Coordinates.geo2xyz(np.radians(M.lon), np.radians(M.lat), M.he)
                # Observation equation in local frame
                
                Ai = np.zeros([3, 3], float)
                Ai[0, 1] = z
                Ai[0, 2] = -y
                Ai[1, 0] = -z
                Ai[1, 2] = x
                Ai[2, 0] = y
                Ai[2, 1] = -x
                
                RAi = np.dot(R, Ai)
                # Predicted velocity
                Pi = np.dot(RAi, W)
                pve = Pi[0, 0] * 1.E3
                pvn = Pi[1, 0] * 1.E3
                
                N = M.copy()
                
                N.Ve = M.Ve - pve
                N.Vn = M.Vn - pvn
                
                lGpoint.append(N)
        
        res_vel = Velocity_Field(lgmt_points=lGpoint)

        return res_vel
        
# !!!!!!!!!!!!! NOT TESTED YET       

###################################################################
    def common(self, linGpoint, prefit=10.0, lexclude=[]):
###################################################################
        """
        Returns a list of sites common to the current Velocity Field object and the list of GMT_Points provided in argument
        Coordinates will be the ones from the Sinex and NOT from the list of Gpoints
        Commons point are points with same code pt and soln
        """
        loutGpoint = []
        lcode = self.lcode()
        for M in linGpoint:
            if M.code not in lcode:continue
            lGpoint = self.subset(lcode=[M.code])
            for N in lGpoint:
                    # Test Prefit coordinates
                if (M.xyz_distance(N) > prefit):
                    print("! Bad prefit (threshold value %8.3lf) for %s: %10.1lf" % (prefit , M.code, M.xyz_distance(N)))
                    break
                else:
                    # 'Match !'
                    loutGpoint.append(N) 
        return(loutGpoint)
    

###################################################################
    def proj_profile(self, slon, slat, elon, elat, d, save=None, verbose=False):
###################################################################
        """
        project velocity components along a great circle defined by initial/stop points (slon,slat,elon,elat)
        
        :param slon,slat: coordinates of profile start point (decimal degrees) 
        :param elon,elat: coordinates of profile end   point (decimal degrees)
        :param d        : maximum distance for a point to be considered
        :param save     : output file name (optional)
        
        :return         :  numpy 1D array with
        np_code, np_distance_along_profile, np_distance_to_profile , \
                np_Ve , np_Vn , np_SVe , np_SVn , \
                np_v_parallele , np_v_perpendicular , \
                np_sigma_v_parallele , np_sigma_v_perpendicular , np_lazimuth
        """

        lcode = []
        ldistance_along_profile = []
        ldistance_to_profile = []

        lVe = []
        lVn = []
        lSVe = []
        lSVn = []
        lv_parallele = []
        lv_perpendicular = []
        lsigma_v_parallele = []
        lsigma_v_perpendicular = []
        lazimuth = []
        
        # Rt in meters
        
        Rt = 6371.0E3
        
        # import
        
        import numpy as np
        from pyacs.lib import coordinates as Coordinates
        from pyacs.lib import vectors
        from pyacs.lib.gmtpoint import GMT_Point
        
        (xi, yi, zi) = Coordinates.geo2xyz(np.radians(slon), np.radians(slat), 0.0)
        Xi = np.array([xi, yi, zi])
        MS = GMT_Point(lon=slon, lat=slat)

        (xs, ys, zs) = Coordinates.geo2xyz(np.radians(elon), np.radians(elat), 0.0)
        Xs = np.array([xs, ys, zs])
        ME = GMT_Point(lon=elon, lat=elat)

        length_arc = MS.spherical_distance(ME) / 1000.0

        POLE = vectors.vector_product(Xi, Xs)
        POLE = vectors.normalize(POLE)

        # LOOP ON SITES

        H_distance = {}
        
        for M in self.sites:
            if verbose:
                print('-- projecting ', M.code)
            (x, y, z) = Coordinates.geo2xyz(np.radians(M.lon), np.radians(M.lat), 0.0)
        
            R = Coordinates.mat_rot_general_to_local(np.radians(M.lon), np.radians(M.lat), unit='radians')
            OM = np.array([x, y, z])
            
            OM = vectors.normalize(OM)
        
            unit_parallele_xyz = vectors.vector_product(POLE, OM)
            unit_parallele_enu = np.dot(R, unit_parallele_xyz)

            v_parallele = M.Ve * unit_parallele_enu[0] + M.Vn * unit_parallele_enu[1]
            v_perpendicular = M.Ve * unit_parallele_enu[1] - M.Vn * unit_parallele_enu[0]
            
            azimuth = np.arctan2(unit_parallele_enu[0], unit_parallele_enu[1])
            
            # longitude of M in the system having AOB as equator and P as north pole
            # to do that, we write the XYZ coordinates of M in the POLE/EQUATOR frame
            # 
            x_m = vectors.scal_prod(OM, vectors.normalize(Xi))
            y_m = vectors.scal_prod(OM, vectors.normalize(vectors.vector_product(Xi, -POLE)))

            longitud_m = np.arctan2(y_m, x_m)
            distance_along_profile = longitud_m * Rt

            # distance to profile is obtained from the latitude in POLE/EQUATOR frame

            z_m = vectors.scal_prod(OM, POLE)
            latitud_m = np.arcsin(z_m)
            
            # distance_to_profile= np.fabs(latitud_m) * Rt
            distance_to_profile = latitud_m * Rt
            
            # now propagates the vcv go get the uncertainty in each direction
            # we use the value of local azimuth of profile and the law of variance propagation
            
            local_R = np.array([[np.cos(azimuth + np.pi / 2.), np.sin(azimuth + np.pi / 2.)], \
                              [-np.sin(azimuth + np.pi / 2.), np.cos(azimuth + np.pi / 2.)]])

            cov = M.SVe * M.SVn * M.SVen
            vcv_en = np.array([[M.SVe ** 2, cov], [cov, M.SVn ** 2]])
            vcv_profile = np.dot(np.dot(local_R, vcv_en), local_R.T)

            sigma_v_parallele = np.sqrt(vcv_profile[0, 0])
            sigma_v_perpendicular = np.sqrt(vcv_profile[1, 1])
            
            if (np.fabs(distance_to_profile / 1000.0) < d) and (distance_along_profile / 1000.0 >= 0.0) and (distance_along_profile / 1000.0 <= length_arc):

                lcode.append(M.code)
                ldistance_along_profile.append(distance_along_profile / 1000.0)
                ldistance_to_profile.append(distance_to_profile / 1000.0)
                lVe.append(M.Ve)
                lVn.append(M.Vn)
                lSVe.append(M.SVe)
                lSVn.append(M.SVn)
                lv_parallele.append(v_parallele)
                lv_perpendicular.append(v_perpendicular)
                lsigma_v_parallele.append(sigma_v_parallele)
                lsigma_v_perpendicular.append(sigma_v_perpendicular)
                lazimuth.append(np.degrees(azimuth))
                
                H_distance[distance_along_profile] = \
                ("%4s           %10.3lf               %10.3lf        %10.2lf %10.2lf %10.2lf %10.2lf    %10.2lf          %10.2lf              %10.2lf            %10.2lf          %10.1lf\n"  \
                % (M.code, distance_along_profile / 1000.0, distance_to_profile / 1000.0, M.Ve, M.Vn , M.SVe, M.SVn, v_parallele, v_perpendicular, sigma_v_parallele, sigma_v_perpendicular, np.degrees(azimuth)))
        
        if save is not None:
            print("-- writing results to ", save)
            fs = open(save, 'w')
            
            fs.write("#site  profile_distance_to_A (km)  distance_to_profile (km)       Ve         Vn        SVe        SVn    V_parallele     V_perpendicular            SV_parallele            SV_perpendicular        azimuth \n")
            fs.write("#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n")
            
            for distance in sorted(H_distance.keys()):
                    fs.write("%s" % H_distance[distance])
                
            fs.close()
                
        np_distance_along_profile = np.array(ldistance_along_profile)
        lindex = np.argsort(np_distance_along_profile)
            
        np_code = np.array(lcode, dtype=str)[lindex]
        np_distance_along_profile = np_distance_along_profile[lindex]
        np_distance_to_profile = np.array(ldistance_to_profile)[lindex]
        np_Ve = np.array(lVe)[lindex]
        np_Vn = np.array(lVn)[lindex]
        np_SVe = np.array(lSVe)[lindex]
        np_SVn = np.array(lSVn)[lindex]
        np_v_parallele = np.array(lv_parallele)[lindex]
        np_v_perpendicular = np.array(lv_perpendicular)[lindex]
        np_sigma_v_parallele = np.array(lsigma_v_parallele)[lindex]
        np_sigma_v_perpendicular = np.array(lsigma_v_perpendicular)[lindex]
        np_lazimuth = np.array(lazimuth)[lindex]

        return  np_code, np_distance_along_profile, np_distance_to_profile , \
                np_Ve , np_Vn , np_SVe , np_SVn , \
                np_v_parallele , np_v_perpendicular , \
                np_sigma_v_parallele , np_sigma_v_perpendicular , np_lazimuth


###################################################################
    def strain(self, lcode, save=None, method='WLS', verbose=False):
###################################################################
        """
        Calculates the strain rate given a list of sites
        
        :param lcode: list of site names 
        :param save  : file to save the result
        :param method: estimator in 'L1','WLS' (default: weighted least-squares 'WLS')
        :param verbose: verbose mode
        
        """

        import numpy as np
        import pyacs.lib.strain
        import pyacs.lib.glinalg
        
        # prepare data
        
        X = np.array([])
        Y = np.array([])
        VE = np.array([])
        VN = np.array([])
        SVE = np.array([])
        SVN = np.array([])
        CVEN = np.array([])
        CODE = np.array([])

        for code in lcode:
            
            M = self.site(code)
            
            X = np.append(X, M.lon)
            Y = np.append(Y, M.lat)
            VE = np.append(VE, M.Ve)
            VN = np.append(VN, M.Vn)
            SVE = np.append(SVE, M.SVe)
            SVN = np.append(SVN, M.SVn)
            CVEN = np.append(CVEN, M.SVen)
            CODE = np.append(CODE, code)
        
        # solve X0,Y0, DV
        
        l0, p0, SOL, COV, RESIDUALS , chi2 = pyacs.lib.strain.vgrad(X, Y, VE, VN, SVE, SVN, CVEN, CODE, method='WLS', verbose=verbose)
        
        CORR , SIGMA = pyacs.lib.glinalg.cov_to_corr(COV)

        [ve0, vn0, dvx_dx, dvx_dy, dvy_dx, dvy_dy] = SOL           
        [sve0, svn0, sdvx_dx, sdvx_dy, sdvy_dx, sdvy_dy] = SIGMA

        sven0 = CORR[0, 1]
        
        # strain rate tensor, rotation rate and shear rates
        
        # From the definition of the strain rate tensor and rotation rate EPS = 0.5 * (DV + DV^T) & omega = 0.5 * (DV - DV^T), 
        # following D. Dong and K. Feigl, we derive EPS, omega and GAMMA using the linear transformation
        
        #     0 0  1  0   0  0  Edot11
        #     0 0  0 .5  .5  0  Edot12
        #     0 0  0  0   0  1  Edot22
        #     0 0  0 .5 -.5  0  omega
        #     0 0  1  0   0 -1  gamma1
        #     0 0  0  1   1  0  gamma2

        TRANS_GRADV = np.zeros((6, 6))
        TRANS_GRADV[0, 2] = 1.
        TRANS_GRADV[1, 3] = .5
        TRANS_GRADV[1, 4] = .5
        TRANS_GRADV[2, 5] = 1.
        TRANS_GRADV[3, 3] = .5
        TRANS_GRADV[3, 4] = -.5
        TRANS_GRADV[4, 2] = 1.
        TRANS_GRADV[4, 5] = -1.
        TRANS_GRADV[5, 3] = 1.
        TRANS_GRADV[5, 4] = 1.
        
        TRANS = np.dot(TRANS_GRADV , SOL)
        
        [edot11, edot12, edot22, omega, gamma1, gamma2] = TRANS
        
        # variance
        
        VAR_TRANS = np.dot(TRANS_GRADV, np.dot(COV , TRANS_GRADV.T))
        
        COV_STRAIN_RATE = VAR_TRANS[:3, :3]
        
        # uncertainties for omega, gamma1 & gamma2
        
        s_omega = np.sqrt(VAR_TRANS[3, 3])
        s_gamma1 = np.sqrt(VAR_TRANS[4, 4])
        s_gamma2 = np.sqrt(VAR_TRANS[5, 5])
        
        # Epsilon1 & 2
        
        eps1 = .5 * (SOL[2] + SOL[5] + np.sqrt((SOL[2] - SOL[5]) ** 2 + (SOL[3] + SOL[4]) ** 2))
        eps2 = .5 * (SOL[2] + SOL[5] - np.sqrt((SOL[2] - SOL[5]) ** 2 + (SOL[3] + SOL[4]) ** 2))
        
        azimut = .5 * np.degrees(np.arctan((SOL[3] + SOL[4]) / (SOL[5] - SOL[2]))) + 90.0
        
        # Gamma
        
        gamma = np.sqrt(gamma1 ** 2 + gamma2 ** 2)
        
        # Uncertainties
        
        term = .5 / np.sqrt((SOL[2] - SOL[5]) ** 2 + (SOL[3] + SOL[4]) ** 2)
        dgdl11 = term * 2. *  (SOL[2] - SOL[5])
        dgdl12 = term * 2. *  (SOL[3] + SOL[4])
        dgdl21 = term * 2. *  (SOL[3] + SOL[4])
        dgdl22 = -term * 2. *  (SOL[2] - SOL[5])
        
        #     partials of epsilon 1 with respect to comps of vel. grad. L
        
        bigg = np.zeros((4, 6))
        #     row 1 corresponds to epsilon 1
        bigg[0, 2] = +1. * (.5 + .5 * dgdl11)
        bigg[0, 3] = +1. * (.5 * dgdl12)
        bigg[0, 4] = +1. * (.5 * dgdl21)
        bigg[0, 5] = +1. * (.5 + .5 * dgdl22)
        
        #     row 2 corresponds to epsilon 2
        bigg[1, 2] = +1. * (.5 - .5 * dgdl11)
        bigg[1, 3] = -1. * (.5 * dgdl12)
        bigg[1, 4] = -1. * (.5 * dgdl21)
        bigg[1, 5] = +1. * (.5 - .5 * dgdl22)
        
        #     row 3 corresponds to theta
        term1 = (SOL[3] + SOL[4]) / (SOL[5] - SOL[2])
        term2 = 1. + term1 ** 2
        term3 = SOL[5] - SOL[2]
        bigg[2, 2] = 0.5 * term1 / term3 / term2
        bigg[2, 3] = 0.5 / term3 / term2
        bigg[2, 4] = 0.5 / term3 / term2
        bigg[2, 5] = -0.5 * term1 / term3 / term2
        
        #     row 4 corresponds to maximum shear
        bigg[3, 2] = dgdl11
        bigg[3, 3] = dgdl12
        bigg[3, 4] = dgdl21
        bigg[3, 5] = dgdl22

        UNC = np.dot(bigg, np.dot(COV, bigg.T))
     
        [seps1, seps2, sazimut, sgamma] = np.sqrt(np.diag(UNC))
        
        # statistics

        dof = 2 * X.size - 6
        
        if dof != 0:
            var_post = np.sqrt(chi2 / dof)
        else:
            var_post = 0.0
        
        wrms = np.sqrt(chi2 / (np.sum(1. / SVE) + np.sum(1. / SVN)))
        
        # significancy
        
        SRV = np.array([edot11, edot12, edot22]).reshape(-1, 1)
        chi2_strain_rate = np.dot(SRV.T, np.dot(pyacs.lib.glinalg.cov_to_invcov(COV_STRAIN_RATE), SRV))
        
        chi2_rotation_rate = omega ** 2 / s_omega ** 2
        
        # print results
        
        RS = ("INFORMATION\n")
        RS += ("----------\n")
        RS += ("file = %s\n" % self.file_name)
        RS += ("Number of sites used in strain rate calculation: %d \n" % (X.size))
        RS += ("%s \n" % (' ').join(lcode))
        RS += ("Velocity at barycenter %10.5lf %10.5lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf\n" % (l0, p0, ve0, vn0, sve0, svn0, sven0))
        
        RS += ("\n")
        
        RS += ("VELOCITY GRADIENT TENSOR VALUES\n")
        RS += ("-------------------------------\n")
        RS += ("  dVx/dx     dVx/dy     dVy/dx      dVy/dy\n")
        RS += ("----------  ---------  ----------  ---------\n")
        RS += ("%10.3lf %10.3lf %10.3lf %10.3lf\n" % (dvx_dx, dvx_dy, dvy_dx, dvy_dy))
        RS += ("VELOCITY GRADIENT TENSOR UNCERTAINTIES\n")
        RS += ("--------------------------------------\n")
        RS += (" sdVx/dx     sdVx/dy    sdVy/dx     sdVy/dy\n")
        RS += ("----------   --------   ---------    --------\n")
        RS += ("%10.3lf %10.3lf %10.3lf %10.3lf\n" % (sdvx_dx, sdvx_dy, sdvy_dx, sdvy_dy))
        RS += ("Unit: nstrain / year \n")
        RS += ("FULL COVARIANCE MATRIX (not rescaled)\n")
        RS += ("-------------------------------------\n")

        RS += ("ve     %10.5lf\n" % (COV[0, 0]))
        RS += ("vn     %10.5lf %10.5lf\n" % tuple(COV[1, :2]))
        RS += ("dvx/dx %10.5lf %10.5lf %10.5lf\n" % tuple(COV[2, :3]))
        RS += ("dvx/dy %10.5lf %10.5lf %10.5lf %10.5lf\n" % tuple(COV[3, :4]))
        RS += ("dvy/dx %10.5lf %10.5lf %10.5lf %10.5lf %10.5lf \n" % tuple(COV[4, :5]))
        RS += ("dvy/dy %10.5lf %10.5lf %10.5lf %10.5lf %10.5lf %10.5lf\n" % tuple(COV[5, :6]))

        RS += ("\n")

        RS += ("STRAIN RATE TENSOR\n")
        RS += ("------------------\n")
        RS += ("     Exx         Exy       Eyy\n")
        RS += ("  --------   --------  ---------\n")
        RS += ("%10.3lf %10.3lf %10.3lf \n" % tuple(TRANS[:3]))

        RS += ("STRAIN RATE COVARIANCE\n")
        RS += ("----------------------\n")
        RS += ("          Exx          Exy           Eyy\n")
        RS += ("Exx %10.5lf \n" % (COV_STRAIN_RATE[0, 0]))
        RS += ("Exy %10.5lf %10.5lf \n" % tuple(COV_STRAIN_RATE[1, :2]))
        RS += ("Eyy %10.5lf %10.5lf %10.5lf \n" % tuple(COV_STRAIN_RATE[2, :3]))
        
        RS += ("\n")
        RS += ("PRINCIPAL AXIS \n")
        RS += ("---------------\n")
        RS += ("\n")
        RS += ("Eps1 = %10.5lf +- %10.5lf \n" % (eps1, seps1))
        RS += ("Eps2 = %10.5lf +- %10.5lf \n" % (eps2, seps2))
        RS += ("Azimut : %10.5lf +- %10.5lf \n" % (azimut, np.degrees(sazimut)))

        RS += ("\n")
        RS += ("Unit: nstrain / year and decimal degree\n")
        RS += ("Eps1: most extensional eigenvalue of strain tensor\n")
        RS += ("Eps2: most compressional eigenvalue of strain tensor\n")
        RS += ("Extension is taken positive\n")
        RS += ("azimut is the one of eps2 in degrees CW from North.\n")

        RS += ("\n")
        RS += ("ROTATION\n")
        RS += ("--------\n")
        RS += ("rotation in deg. / Myr :  %10.3lf +- %10.3lf \n" % (omega * 1.E-9 * 180. / np.pi * 1.E6, s_omega * 1.E-9 * 180. / np.pi * 1.E6))
        RS += ("rotation in rad /  yr  :  %.4E +- %.4E \n" % (omega * 1.E-9, s_omega * 1.E-9))

        RS += ("\n")
        RS += ("SHEAR RATES\n")
        RS += ("-----------\n")
        RS += ("gamma1 in nstrain / yr : %10.5lf +- %10.5lf \n" % (gamma1, s_gamma1))
        RS += ("gamma2 in nstrain / yr : %10.5lf +- %10.5lf \n" % (gamma2, s_gamma2))
        RS += ("gamma  in nstrain / yr : %10.5lf +- %10.5lf \n" % (gamma, sgamma))

        RS += ("\n")
        RS += ("STATISTICS\n")
        RS += ("----------\n")
        RS += ("chi2                       : %10.5lf \n" % (chi2))
        RS += ("degree of freedom          : %10d \n" % (dof))
        RS += ("posterior variance factor  : %10.5lf \n" % (var_post))
        RS += ("wrms (mm/yr)               : %10.5lf \n" % (wrms))

        RS += ("\n")
        RS += ("TEST OF SIGNIFICANCY\n")
        RS += ("--------------------\n")
        RS += ("\n")
        
        RS += ("tested              chi2      threshold values\n")
        RS += ("parameter           values       95%       99%\n")
        RS += ("-----------------------------------------------\n")
        RS += ("strain rate tensor %10.2lf   7.82     11.35\n" % chi2_strain_rate)
        RS += ("rotation rate      %10.2lf   3.84      6.64\n" % chi2_rotation_rate)

        RS += ("\n")
        RS += ("RESIDUALS\n")
        RS += ("---------\n")
        RS += ("code      long.       lat.       Ve        Vn    Res_Ve    Res_Vn   sigma_Ve   sigma_Vn  RN_Ve     RN_Vn\n")
        RS += ("---------------------------------------------------------------------------------------------------------\n")
       
        for i in np.arange(X.size):
            RS += ("%s %10.5lf %10.5lf %8.2lf  %8.2lf  %8.2lf  %8.2lf  %8.2lf  %8.2lf  %8.2lf  %8.2lf  \n" % \
                   (CODE[i] , X[i] , Y[i], VE[i], VN[i], RESIDUALS[2 * i] , RESIDUALS[2 * i + 1] , SVE[i], SVN[i] , RESIDUALS[2 * i] / SVE[i], RESIDUALS[2 * i + 1] / SVN[i]))
        
        if verbose or (save is None):
            print(RS)
        
        if save is not None:
            f = open(save, 'w')
            f.write(RS)
            f.close()
        
        return self
        
