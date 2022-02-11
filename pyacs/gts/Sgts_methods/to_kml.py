def to_kml( self , kml_file ):
    """
    write a kml file with site location

    :param kml_file: output kml file
    :return: self
    """
    # adapted from https://gist.github.com/mazzma12/0a32ce693bb42b742252caabb98519db

    import fiona
    import geopandas as gpd
    from shapely.geometry import Point

    # Enable fiona driver
    gpd.io.file.fiona.drvsupport.supported_drivers['KML'] = 'rw'

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
    d = {'NAME': lcode, 'geometry': lPoint}
    gdf = gpd.GeoDataFrame(d, crs="EPSG:4326")

    # Write file
    with fiona.Env():
        # Might throw a WARNING - CPLE_NotSupported in b'dataset sample_out.kml does not support layer creation option ENCODING'
        gdf.to_file( kml_file, driver='KML')
