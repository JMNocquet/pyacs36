"""Convert a pyeblock fault file into a shapefile."""


def _add_field(w, name, field_type, size, decimal=0):
    """Add a shapefile field; compatible with pyshp 2.x (string types) and 3.x (no 'I' alias)."""
    # pyshp 3.x: no 'I' (integer), use 'N' with decimal=0; expects (name, type, size, decimal).
    if field_type == 'I':
        field_type = 'N'
    try:
        w.field(name, field_type, size, decimal)
    except TypeError:
        # pyshp 2.x may only accept 3 args (name, type, size)
        w.field(name, field_type, size)


def pyeblock_fault(pyeblock_file, shp_name, verbose=False):
    """Convert a pyeblock fault file into a shapefile.

    Parameters
    ----------
    pyeblock_file : str
        Path to the pyeblock fault file.
    shp_name : str
        Output shapefile base name (no extension).
    verbose : bool, optional
        If True, print progress. Default is False.
    """
    import shapefile
    import numpy as np

    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR

    lfaults = []
    lrecord = []

    try:
        FAULTS = np.array(np.asmatrix(np.genfromtxt(pyeblock_file, usecols=range(8), comments='#')))
        FAULTS_POLY = np.array(np.asmatrix(np.genfromtxt(pyeblock_file, usecols=(8, 9), dtype=str)))
    except Exception:
        ERROR("Could not read %s" % pyeblock_file)

    w = shapefile.Writer(shp_name, shapeType=shapefile.POLYLINE)
    # pyshp 3.x removed 'I' (integer) alias; use 'N' (numeric) with decimal=0
    _add_field(w, 'ID', 'N', 40, 0)
    _add_field(w, 'top', 'F', 40, 6)
    _add_field(w, 'bottom', 'F', 40, 6)
    _add_field(w, 'dip', 'F', 40, 6)
    _add_field(w, 'left_block', 'C', 40, 0)
    _add_field(w, 'right_block', 'C', 40, 0)

    for i in np.arange(FAULTS.shape[0]):
        VERBOSE(str(FAULTS[i, 1:3]))
        w.line([[FAULTS[i, 1:3], FAULTS[i, 3:5]]])
        w.record(str(i), FAULTS[i, 5], FAULTS[i, 6], FAULTS[i, 7], FAULTS_POLY[i, 0], FAULTS_POLY[i, 1])

    shp_path = shp_name + '.shp'
    VERBOSE(("saving shapefile %s " % (shp_path)))
    w.close()

    prj = open("%s.prj" % shp_name, "w")
    epsg = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]'
    prj.write(epsg)
    prj.close()

    gmtfile = shp_name + '.gmt'
    VERBOSE(("saving gmt file %s " % gmtfile.split('/')[-1]))
    f = open(gmtfile, 'w')
    for i in np.arange(FAULTS.shape[0]):
        f.write('> -Z%d\n' % i)
        f.write("%10.3lf %10.3lf\n" % (FAULTS[i, 1], FAULTS[i, 2]))
        f.write("%10.3lf %10.3lf\n" % (FAULTS[i, 3], FAULTS[i, 4]))
    f.write('>\n')
    f.close()
