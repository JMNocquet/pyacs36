"""Convert a GMT psvelo file into a shapefile."""


def psvelo_to_shapefile(psvelo_file, shp_name, verbose=False):
    """Convert a GMT psvelo file into a shapefile.

    Parameters
    ----------
    psvelo_file : str
        Path to the GMT psvelo file.
    shp_name : str
        Output shapefile base name (no extension).
    verbose : bool, optional
        If True, print progress. Default is False.
    """
    import shapefile
    import numpy as np
    from tqdm import tqdm

    import logging
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.debug_message as DEBUG

    VERBOSE("Reading %s" % psvelo_file)

    try:
        np_vel = np.atleast_2d(np.genfromtxt(psvelo_file, comments='#', usecols=(0, 1, 2, 3, 4, 5, 6)))
        np_vel_name = np.atleast_1d(np.genfromtxt(psvelo_file, comments='#', usecols=(7), dtype=str))
    except Exception:
        ERROR("Could not read %s" % psvelo_file, exit=True)

    if np_vel.size == 0:
        return

    VERBOSE("%d records read" % np_vel.shape[0])

    VERBOSE("Initialize shapefile")
    w = shapefile.Writer(shp_name, shapeType=shapefile.POINT)
    w.field('longitude', 'C', '40')
    w.field('latitude', 'C', '40')
    w.field('Ve', 'C', '40')
    w.field('Vn', 'C', '40')
    w.field('S_Ve', 'C', '40')
    w.field('S_Vn', 'C', '40')
    w.field('S_Ven', 'C', '40')
    w.field('name', 'C', '40')
    w.field('name_velocity', 'C', '40')

    VERBOSE("Populating shapefile")

    if logging.getLogger("my_logger").getEffectiveLevel() in [logging.INFO, logging.DEBUG]:
        for i in tqdm(np.arange(np_vel.shape[0]), desc='Filling shapefile'):
            DEBUG(("adding %s " % np_vel_name[i]))
            w.point(np_vel[i, 0], np_vel[i, 1])
            name_velocity = ("%s %5.2lf" % (np_vel_name[i], np.sqrt(np_vel[i, 2]**2 + np_vel[i, 3]**2)))
            w.record(np_vel[i, 0], np_vel[i, 1], np_vel[i, 2], np_vel[i, 3], np_vel[i, 4], np_vel[i, 5], np_vel[i, 6], np_vel_name[i], name_velocity)
    else:
        for i in np.arange(np_vel.shape[0]):
            DEBUG(("adding %s " % np_vel_name[i]))
            w.point(np_vel[i, 0], np_vel[i, 1])
            name_velocity = ("%s %5.2lf" % (np_vel_name[i], np.sqrt(np_vel[i, 2]**2 + np_vel[i, 3]**2)))
            w.record(np_vel[i, 0], np_vel[i, 1], np_vel[i, 2], np_vel[i, 3], np_vel[i, 4], np_vel[i, 5],
                    np_vel[i, 6], np_vel_name[i], name_velocity)

    shp_path = shp_name + '.shp'
    w.close()
    VERBOSE("saving shapefile %s " % (shp_path))

    DEBUG("writing EPSG")
    prj = open("%s.prj" % shp_name, "w")
    epsg = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]'
    prj.write(epsg)
    prj.close()
