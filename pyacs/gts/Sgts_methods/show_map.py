def show_map( self , bounds = None , highlight=[] , geotiff = None, tile = False, grid=True ,show=True, save=False ):
    """

    :param self: Sgts instance
    :param bounds: map bounds as list [min_lon,max_lon,min_lat,max_lat]
    :param highlight: list of site code to be highlighted
    :param geotiff: a global [-180,180,-90,90] geotiff file
    :param tile: boolean. If True reads tiles from Stamen Design using contextily. Requires internet connection. Default is False.
    :param grid: boolean. plot grid and axis labels (default grid=True)
    :return: self
    """
    import geopandas
    import matplotlib.pyplot as plt
    from matplotlib.image import imread
    import numpy as np
    import contextily as ctx
    import cartopy.crs as ccrs
    #import georaster


    # get coordinates
    llong = []
    llat  = []
    lcode = []
    for code in sorted( self.lcode() ):
        lcode.append( code )
        llong.append( self.__dict__[code].lon )
        llat.append( self.__dict__[code].lat )



    # get bounds if not provided
    if bounds is None:
        lon_min = np.min( llong )
        lon_max = np.max( llong )
        lat_min = np.min( llat )
        lat_max = np.max( llat )
        delta_lon = (lon_max-lon_min) / 20.
        delta_lat = (lat_max-lat_min) / 20.
        lon_min = lon_min - delta_lon
        lon_max = lon_max + delta_lon
        lat_min = lat_min - delta_lat
        lat_max = lat_max + delta_lat

    else:
        [lon_min,lon_max,lat_min,lat_max] = bounds
        llong.append( lon_min )
        llong.append( lon_max )
        llat.append( lat_min )
        llat.append( lat_max )

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

    # background image
    if geotiff is not None:
        img = imread(geotiff)
        ax.imshow(img, origin='upper', extent=[-180, 180, -90, 90], transform=ccrs.PlateCarree())
    else:
        # background tile
        if tile:
            ax.plot(llong, llat, 'ro', markersize=3)
            ctx.add_basemap(ax, crs='EPSG:4326')
        else:
            world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
            world.plot( ax=ax )


    ax.set_xlim(lon_min,lon_max)
    ax.set_ylim(lat_min,lat_max)
    ax.plot(llong,llat,'ro',markersize=3)

    for i, txt in enumerate(lcode):
        ax.annotate(txt, (llong[i], llat[i]), fontsize=6)

    for code in highlight:
        ax.annotate(code, (self.__dict__[code].lon, self.__dict__[code].lat), fontsize=8 , color='r')
        ax.plot( [self.__dict__[code].lon], [self.__dict__[code].lat], 'ws',markersize=5)


    # grid
    if grid:
        ax.gridlines(draw_labels=True, dms=False, x_inline=False, y_inline=False)

    # show
    if show:
        plt.ion()
        plt.show()
    else:
        plt.ioff()

    # save
    if save:
        fig.savefig( save )

    return self
