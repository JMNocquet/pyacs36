def show_map( self , bounds=None, highlight=[], tile=False, tile_provider='OpenTopoMap',
              figsize=(7, 9), show=True,
              save=False):
    """

    :param self: Sgts instance
    :param bounds: map bounds as list [min_lon,max_lon,min_lat,max_lat]
    :param highlight: list of site code to be highlighted
    :param tile: boolean. If True reads tiles from tile_provider. Requires internet connection. Default is False.
    :param tile_provider: 'OpenTopoMap' or 'OpenStreetMap.Mapnik'
    :return: self
    """

    import os
    import geopandas
    import matplotlib.pyplot as plt
    from matplotlib.image import imread
    from shapely.geometry import Point
    import numpy as np
    import contextily as ctx
    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG

    plt.ioff()

    # highlight
    if isinstance(highlight, str):
        highlight = [highlight]

    if highlight != []:
        hts = self.sub(linclude=highlight)

    # get coordinates from ts
    llong = []
    llat = []
    lcode = []
    lPoint = []

    for code in sorted(self.lcode()):
        lcode.append(code)
        llong.append(self.__dict__[code].lon)
        llat.append(self.__dict__[code].lat)
        lPoint.append(Point(self.__dict__[code].lon, self.__dict__[code].lat))

    # get bounds if not provided
    if bounds is None:
        lon_min = np.min(llong)
        lon_max = np.max(llong)
        lat_min = np.min(llat)
        lat_max = np.max(llat)
        delta_lon = (lon_max - lon_min) / 20.
        delta_lat = (lat_max - lat_min) / 20.
        lon_min = lon_min - delta_lon
        lon_max = lon_max + delta_lon
        lat_min = lat_min - delta_lat
        lat_max = lat_max + delta_lat

    else:
        [lon_min, lon_max, lat_min, lat_max] = bounds
        llong.append(lon_min)
        llong.append(lon_max)
        llat.append(lat_min)
        llat.append(lat_max)

    d = {'col1': lcode, 'geometry': lPoint}
    gdf = geopandas.GeoDataFrame(d, crs="EPSG:4326")

    ### TILE CASE ###
    if tile:
        wmgdf = gdf.to_crs(3857)
        ax = wmgdf.plot(markersize=5, color='r', figsize=figsize)

        VERBOSE("Loading tiles")
        try:
            if tile_provider == 'OpenTopoMap':
                ctx.add_basemap(ax=ax, source=ctx.providers.OpenTopoMap, alpha=0.5)
            if tile_provider == 'OpenStreetMap.Mapnik':
                ctx.add_basemap(ax=ax, source=ctx.providers.OpenStreetMap.Mapnik)

        except:
            WARNING('Some errors were raised during tiles handling, most probably due to internet connection')
            MESSAGE('Continuing without tiles')
            plt.close()
            tile = False
            pass

    if tile:
        ax.set_title("%d sites" % self.n())
        for i, txt in enumerate(lcode):
            ax.annotate(txt, (wmgdf['geometry'][i].x, wmgdf['geometry'][i].y), fontsize=8)
        ax.axis('off')

        if highlight != []:
            # get coordinates from ts
            llong = []
            llat = []
            lcode = []
            lPoint = []
            for code in sorted(hts.lcode()):
                lcode.append(code)
                llong.append(self.__dict__[code].lon)
                llat.append(self.__dict__[code].lat)
                lPoint.append(Point(self.__dict__[code].lon, self.__dict__[code].lat))
            d = {'col1': lcode, 'geometry': lPoint}
            gdf = geopandas.GeoDataFrame(d, crs="EPSG:4326")
            wmgdf = gdf.to_crs(3857)
            wmgdf.plot(ax=ax, markersize=500, color='b', marker='*', zorder=10)

    if not tile:

        world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
        ax = world.plot(figsize=figsize, alpha=0.3)

        # plot GPS as dots
        ax = gdf.plot(ax=ax, markersize=5, color='r')

        ax.set_xlim(lon_min, lon_max)
        ax.set_ylim(lat_min, lat_max)
        # ax.plot(llong,llat,'ro',markersize=3)

        ax.set_title("%d sites" % self.n())

        for i, txt in enumerate(lcode):
            ax.annotate(txt, (llong[i], llat[i]), fontsize=8)

        for code in highlight:
            if code in self.lcode():
                ax.annotate(code, (self.__dict__[code].lon, self.__dict__[code].lat), fontsize=8, color='r')
                ax.plot([self.__dict__[code].lon], [self.__dict__[code].lat], \
                        marker='*', color='yellow', markersize=20, markeredgecolor='k', zorder=10)

        # grid
        ax.grid(color='grey')

        # draw_labels=True, dms=False, x_inline=False, y_inline=False

    # show
    if show:
        plt.ion()
        plt.show()
    else:
        plt.ioff()

    # save
    if save:
        plt.savefig(save)
        VERBOSE("Saved plot as %s " % (os.path.abspath(save)))

    return self