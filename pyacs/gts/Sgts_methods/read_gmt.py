###################################################################
def read_gmt(self,gmt=True,verbose=False,vel=False):
###################################################################
    """
    Reads a gmt psvelo file and populates lon and lat attributes for each Gts of Sgts
     
    :param gmt: if True tries to read '/../stat/pyacs_void.gmt', if a string then it is the gmt file to be read. 
    :param verbose: verbose mode
    :param vel: boolean. If True, fills the .velocity attribute of every time series with the values read in the gmt file.
     
    :note: this method is always in place.
    """
    
    if gmt==True:
        gmt_file=self.dir+'/../stat/pyacs_void.gmt'
    if  isinstance(gmt,str):
        gmt_file=gmt
     
    if gmt != False:
        if verbose:print("-- Importing modules Velocity_Field, Gpoint ",gmt_file)
    
        from pyacs.lib.vel_field import Velocity_Field
         
        if verbose:print("-- Reading ",gmt_file)
        velf=Velocity_Field()
        velf.read(file_name=gmt_file)
        lsite=vel.lcode()
         
        for gts in self.lGts():
            if verbose:print("-- Processing ",gts.code)
            if gts.code in lsite:
                M=velf.site(gts.code)
                gts.lon=M.lon
                gts.lat=M.lat
                if vel:
                    gts.velocity=[M.Ve,M.Vn,0.0,M.SVe,M.SVn,0.0]
