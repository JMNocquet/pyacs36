def show_map(self, bounds=None, resolution='l', highlight=[], tile=False, tile_provider='OpenTopoMap',
             figsize=(7, 9), show=True,
             save=False, velocity=False, title=None):
    """Show map of Sgts sites (optionally with velocity).

    Parameters
    ----------
    bounds : list, optional
        [min_lon, max_lon, min_lat, max_lat]. Default is None (auto).
    resolution : str, optional
        'c', 'l', 'i', 'h', 'f'. Default is 'l'.
    highlight : list, optional
        Site codes to highlight. Default is [].
    tile : bool, optional
        If True, use tile_provider (requires internet). Default is False.
    tile_provider : str, optional
        'OpenTopoMap' or 'OpenStreetMap.Mapnik'. Default is 'OpenTopoMap'.
    figsize : tuple, optional
        Figure size. Default is (7, 9).
    show : bool, optional
        If True, show plot. Default is True.
    save : bool, optional
        If True, save figure. Default is False.
    velocity : bool, optional
        If True, plot velocities. Default is False.
    title : str, optional
        Plot title. Default is None.

    Returns
    -------
    self
    """

    import os
    import geopandas
    import matplotlib.pyplot as plt
    from matplotlib.image import imread
    from shapely.geometry import Point
    import numpy as np
    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG

    import inspect
    VERBOSE("Running Sgts.%s" % inspect.currentframe().f_code.co_name)

    # handle error
    if self.n() == 0:
        ERROR("Sgts has no time series")

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
        # check if contextily is installed
        try:
            import contextily as ctx
        except:
            ERROR("tile option requires contextily which is not installed. Please install it with 'pip install contextily'")
            return self
        else:
            pass
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
        ax.set_title(title + " - %d sites" % self.n() if title is not None else "%d sites" % self.n())
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

        if velocity:
            np_arrow = np.zeros((self.n(), 4))
            for i, code in enumerate(self.lcode()):
                vel = self.__dict__[code].velocity * 1.E3
                np_arrow[i, 0] = self.__dict__[code].lon
                np_arrow[i, 1] = self.__dict__[code].lat
                np_arrow[i, 2] = vel[1]
                np_arrow[i, 3] = vel[0]

            ax.quiver(np_arrow[:, 0], np_arrow[:, 1], np_arrow[:, 2], np_arrow[:, 3])

    ### NO TILE CASE ###
    if not tile:
        # deprecated since geopandas stopped providing the world map
        # world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
        # ax = world.plot(figsize=figsize, alpha=0.3)

        dir_data = os.path.expanduser('~') + '/.pyacs/data/GMT'
        wshp = dir_data + '/GSHHS_shp/' + resolution + '/GSHHS_' + resolution + '_L1.shp'

        if not os.path.exists(wshp):

            # creates $HOME/.pyacs if it does not exists
            from pathlib import Path
            Path(dir_data).mkdir(parents=True, exist_ok=True)

            url = 'http://www.soest.hawaii.edu/pwessel/gshhg/gshhg-shp-2.3.7.zip'
            VERBOSE("Missing shapefile for background map")
            VERBOSE("Downloading world map: %s" % url)
            VERBOSE("World map will be stored in: %s" % os.path.expanduser('~') + '/.pyacs/data/GMT')

            import requests
            from tqdm import tqdm
            import zipfile

            def download(url: str, fname: str):
                resp = requests.get(url, stream=True)
                total = int(resp.headers.get('content-length', 0))
                # Can also replace 'file' with a io.BytesIO object
                with open(fname, 'wb') as file, tqdm(
                        desc=fname,
                        total=total,
                        unit='iB',
                        unit_scale=True,
                        unit_divisor=1024,
                ) as bar:
                    for data in resp.iter_content(chunk_size=1024):
                        size = file.write(data)
                        bar.update(size)

            def unzip_file(zip_filepath, dest_dir):
                with zipfile.ZipFile(zip_filepath, 'r') as zip_ref:
                    zip_ref.extractall(dest_dir)

            try:
                download(url=url, fname=dir_data + '/gshhg-shp-2.3.7.zip')
            except:
                ERROR("Could not download %s" % url)

            try:
                unzip_file(zip_filepath=dir_data + '/gshhg-shp-2.3.7.zip', dest_dir=dir_data)
            except:
                ERROR("Could not unzip %s" % dir_data + '/gshhg-shp-2.3.7.zip')

            # clean
            try:
                os.remove(dir_data + '/gshhg-shp-2.3.7.zip')
            except:
                ERROR("Could not remove %s" % dir_data + '/gshhg-shp-2.3.7.zip')

        # read shapefile
        VERBOSE("Reading %s" % wshp)
        world = geopandas.read_file(wshp)
        ax = world.plot(figsize=figsize, alpha=0.3)
        # plot GPS as dots
        ax = gdf.plot(ax=ax, markersize=5, color='r')

        ax.set_xlim(lon_min, lon_max)
        ax.set_ylim(lat_min, lat_max)
        # ax.plot(llong,llat,'ro',markersize=3)

        ax.set_title(title if title is not None else "%d sites" % self.n())

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

        if velocity:
            np_arrow = np.zeros((self.n(), 4))
            for i, code in enumerate(self.lcode()):
                vel = self.__dict__[code].velocity * 1.E3
                np_arrow[i, 0] = self.__dict__[code].lon
                np_arrow[i, 1] = self.__dict__[code].lat
                np_arrow[i, 2] = vel[1]
                np_arrow[i, 3] = vel[0]

            ax.quiver(np_arrow[:, 0], np_arrow[:, 1], np_arrow[:, 2], np_arrow[:, 3])

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