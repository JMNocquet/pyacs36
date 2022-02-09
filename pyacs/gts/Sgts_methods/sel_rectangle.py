###################################################################
def sel_rectangle(self,bounds,verbose=True):
###################################################################
    """
    selects the time series for sites within a rectangles
     
    :param bounds: [lon_min,lon_max,lat_min,lat_max]
    :pram verbose: verbose mode
     
    :return: a new Sgts instance
    """


    # import    
    from pyacs.gts.Sgts import Sgts
    
    
    [lon_min,lon_max,lat_min,lat_max]=bounds
    
    new_Sgts=Sgts(read=False)
    for gts in self.lGts():
        current_lon=gts.lon
        current_lat=gts.lat
         
        if current_lon>=lon_min and current_lon<=lon_max<=lon_max and current_lat>= lat_min and current_lat<=lat_max:
            if verbose:print("-- ",gts.code," selected")
            new_Sgts.append(gts)

    return(new_Sgts)
