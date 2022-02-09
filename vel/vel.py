class vel:
    
    """
    vel class: horizontal velocity field class
    
    A velocity field is made of pyacs.lib.gmtpoint objects
    
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
