"""Grid generation utilities."""


def make_grid(min_lon, max_lon, min_lat, max_lat,
              nx=None, ny=None, step_x=None, step_y=None, outfile=None, format='psvelo', comment=''):
    """Generate a regular grid and optionally write it to a text file.

    Parameters
    ----------
    min_lon : float
        Minimum longitude (decimal degrees).
    max_lon : float
        Maximum longitude (decimal degrees).
    min_lat : float
        Minimum latitude (decimal degrees).
    max_lat : float
        Maximum latitude (decimal degrees).
    nx : int, optional
        Number of points along longitude. If provided, ny defaults to nx.
    ny : int, optional
        Number of points along latitude.
    step_x : float, optional
        Longitude step (alternative to nx).
    step_y : float, optional
        Latitude step; defaults to step_x if not provided.
    outfile : str, optional
        Output file path. If None, no file is written.
    format : str, optional
        None for lon,lat only; 'psvelo' for psvelo-style rows. Default is 'psvelo'.
    comment : str, optional
        Comment for the output file. Default is ''.

    Returns
    -------
    numpy.ndarray
        2D array of shape (n_points, 2) with lon, lat columns.
    """
    import numpy as np
    from colors import red

    if [nx, step_x] == [None, None]:
        error_str = red("[PYACS ERROR] please provide a value for argument nx or step_x")
        print(error_str)
        raise TypeError

    if (nx is not None) and (step_x is not None):
        error_str = red("[PYACS ERROR] nx or step_x argument need to be specified, not both")
        print(error_str)
        raise TypeError

    if nx is not None:
        if ny is None:
            ny = nx
        x = np.linspace(min_lon, max_lon, nx)
        y = np.linspace(min_lat, max_lat, ny)

    if step_x is not None:
        if step_y is None:
            step_y = step_x
        x = np.arange(min_lon, max_lon, step_x)
        y = np.arange(min_lat, max_lat, step_y)

    mesh_grid = np.meshgrid(x, y)
    grid = np.vstack((mesh_grid[0].T.flatten(), mesh_grid[1].T.flatten())).T

    if format == 'psvelo':
        output_np = np.zeros((grid.shape[0], 8))
        output_np[:, :2] = grid
        output_np[:, 7] = np.arange(grid.shape[0])

        if outfile is not None:
            np.savetxt(outfile, output_np, "%10.6lf %10.6lf %4.1lf %4.1lf %4.1lf %4.1lf %4.1lf %04d", header=comment)

    else:
        if outfile is not None:
            np.savetxt(outfile, grid, "%10.6lf %10.6lf", header=comment)

    return grid

