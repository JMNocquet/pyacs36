"""
Routine to force daily solution at 12:00:00
"""

###################################################################
def force_daily( self , in_place=False ):
###################################################################
    """
    force a time series to be daily with dates at 12:00:00 of every day
    """

    import pyacs.lib.astrotime

    # decyear -> mjd -> int + 0.5    
    np_coor_mjd = pyacs.lib.astrotime.decyear2mjd(self.data[:,0]).astype(int) + 0.5
    np_coor_decyear = pyacs.lib.astrotime.mjd2decyear(np_coor_mjd)

    if in_place:
        self.data[:,0] = np_coor_decyear
        if self.data_xyz is not None:
            self.data_xyz[:,0] = np_coor_decyear
            
        return(self)
    else:
        new_Gts=self.copy()

        new_Gts.data[:,0] = np_coor_decyear
        
        if new_Gts.data_xyz is not None:
            new_Gts.data_xyz[:,0] = np_coor_decyear
        return(new_Gts)

    

