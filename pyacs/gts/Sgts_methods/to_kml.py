def to_kml( self , kml_file ):
    """Write a KML file with site locations.

    Parameters
    ----------
    kml_file : str
        Output KML file path.

    Returns
    -------
    self
    """
    # adapted from https://gist.github.com/mazzma12/0a32ce693bb42b742252caabb98519db

    import fiona
    import geopandas as gpd
    from shapely.geometry import Point

    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG

    import inspect
    VERBOSE("Running Sgts.%s" % inspect.currentframe().f_code.co_name)

    def write_kml_from_coords(lon, lat, h, name, output_path="output.kml"):
        """
        Writes a simple KML file with placemarks from coordinate and name lists.

        Parameters:
            lon (list): List of longitudes
            lat (list): List of latitudes
            h (list): List of altitudes (heights)
            name (list): List of placemark names
            output_path (str): Output KML filename
        """
        if not (len(lon) == len(lat) == len(h) == len(name)):
            raise ValueError("All input lists must be the same length")

        kml_header = """<?xml version="1.0" encoding="UTF-8"?>
    <kml xmlns="http://www.opengis.net/kml/2.2">
    <Document>
    """
        kml_footer = "</Document>\n</kml>\n"

        placemarks = ""
        for i in range(len(lon)):
            placemarks += f"""  <Placemark>
        <name>{name[i]}</name>
        <Point>
          <coordinates>{lon[i]},{lat[i]},{h[i]}</coordinates>
        </Point>
      </Placemark>
    """

        kml_content = kml_header + placemarks + kml_footer

        with open(output_path, "w", encoding="utf-8") as f:
            f.write(kml_content)

        VERBOSE(f"KML file written to: {output_path}")

    # get coordinates from ts
    llong = []
    llat = []
    lh = []
    lcode = []

    for code in sorted(self.lcode()):
        lcode.append(code)
        llong.append(self.__dict__[code].lon)
        llat.append(self.__dict__[code].lat)
        lh.append(self.__dict__[code].h)

    write_kml_from_coords(llong, llat, lh, lcode, kml_file)

    return self
