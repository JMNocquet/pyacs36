def show_ivel_map_gmt(self, day,
                        topo=False,
                        projection="M10",
                        region=None,
                        eq_lon_lat=None,
                        scale_gps=0.1,
                        gps_scale=0.01,
                        no_clip=True,
                        show='external',
                        save_dir=None):
    """
    Make a map of instantaneous velocity (ivel) using pygmt.

    Parameters:
    day (tuple, datetime, float, or list): Day to read. Can be (day, month, year) tuple, datetime object, decimal year, or a list of these.
    topo (bool): If True, shaded topography will be shown in the background. Default is False.
    projection (str): Projection type for the map. Default is "M10".
    region (list): List defining the region to plot [min_lon, max_lon, min_lat, max_lat]. If None, it will be determined from data.
    eq_lon_lat (tuple): Tuple of earthquake longitude and latitude to plot. Default is None.
    scale_gps (float): Scale for GPS vectors. Default is 0.1.
    gps_scale (float): Scale for GPS uncertainty vectors. Default is 0.01.
    no_clip (bool): If True, vectors will not be clipped at the map boundaries. Default is True.
    show (str): 'external' or 'notebook' to show the plot. Anything else will not show the plot. Default is 'external'.
    save_dir (str): Directory to save the plot. If None, the plot will not be saved. Default is None.

    Returns:
    None
    """
    import numpy as np
    import pygmt
    import pyacs.lib.astrotime as at
    from datetime import datetime

    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG
    import pyacs.debug
    import pandas as pd
    from pyacs.gts.lib.tensor_ts.sgts2obs_tensor import sgts2tensor

    def decipher_day(day):
        """
        Decipher the input day and convert it into integer MJD.

        Parameters:
        day (tuple, datetime, float): Day to decipher.

        Returns:
        int: MJD of the input day.
        """
        try:
            if isinstance(day, tuple):
                return int(at.cal2mjd(day[0], day[1], day[2]))
            elif isinstance(day, float) and day > 1980:
                return int(at.decyear2mjd(day))
            elif isinstance(day, datetime):
                return int(at.datetime2mjd(day))
            else:
                raise ValueError("Invalid day format")
        except Exception as e:
            ERROR(f"Could not decipher day: {day}. Error: {e}", exit=True)

    lmjd = [decipher_day(d) for d in day] if isinstance(day, list) else [decipher_day(day)]

    # Convert Sgts to obs_tensor
    OBS, np_names, np_coor, np_seconds_sol = sgts2tensor(self)
    np_coor = np_coor[:, :2]
    np_mjd = np.array(at.datetime2mjd(at.seconds2datetime(np_seconds_sol)), dtype=int)

    # Map settings and pygmt parameters
    arrow_shape = ".4c+n50/0.5+h.5+a40+gblack+e+z+qk"

    if topo:
        grid = pygmt.datasets.load_earth_relief(resolution="01m", region=region)
        shade = pygmt.grdgradient(grid=grid, azimuth="0/90", normalize="t1")

    if region is None:
        min_lon, min_lat = np.min(np_coor, axis=0)
        max_lon, max_lat = np.max(np_coor, axis=0)
        delta_lon = max_lon - min_lon
        delta_lat = max_lat - min_lat
        region = [min_lon - 0.05 * delta_lon, max_lon + 0.05 * delta_lon, min_lat - 0.05 * delta_lat, max_lat + 0.05 * delta_lat]

    arrow_scale = f"e{gps_scale:.3f}/0./0"

    TMP_OBS = np.copy(OBS)
    TMP_OBS[np.isnan(TMP_OBS)] = 0

    for mjd_ivel in lmjd:
        str_date = at.mjd2datetime(mjd_ivel).isoformat()[:10]

        try:
            idx = np.where(np_mjd == mjd_ivel)[0][0]
        except IndexError:
            ERROR(f"No day {str_date} in the current Sgts", exit=False)
            return

        lidx_largest = np.argsort(TMP_OBS[idx, :, 0]**2 + TMP_OBS[idx, :, 1]**2)[-5:]
        print("-------------------------------------------------------")
        for i in lidx_largest:
            print(f"{str_date} {np_names[i]} {OBS[idx, i, 0]:10.1f} {OBS[idx, i, 1]:10.1f}")

        df = pd.DataFrame({
            "x": np_coor[:, 0],
            "y": np_coor[:, 1],
            "east_velocity": OBS[idx, :, 0],
            "north_velocity": OBS[idx, :, 1],
            "east_sigma": np.ones(np_names.shape[0]),
            "north_sigma": np.ones(np_names.shape[0]),
            "correlation_EN": np.zeros(np_names.shape[0]),
            "SITE": np_names
        })

        df_scale = pd.DataFrame({
            "x": region[0] + 0.3 * (region[1] - region[0]),
            "y": region[3] - 0.02 * (region[3] - region[2]),
            "east_velocity": [100],
            "north_velocity": [0],
            "east_sigma": [1],
            "north_sigma": [1],
            "correlation_EN": [1],
            "SITE": '100 mm/yr'
        })

        fig = pygmt.Figure()

        pygmt.makecpt(cmap="terra")
        pygmt.config(MAP_FRAME_TYPE="plain")
        pygmt.config(FORMAT_GEO_MAP="ddd.x")
        pygmt.config(FONT_TITLE=14)

        title = f"+tHorizontal daily velocity {str_date}"
        fig.basemap(region=region, projection=projection, frame=["a1f1", title])
        if topo:
            fig.grdimage(grid=grid, cmap=True, shading=shade, transparency=50)
            fig.coast(shorelines=1)
        else:
            fig.coast(land="grey80", water="191/239/255")

        # check that df.east_velocity and df.north_velocity are not all nan
        if np.isnan(df.east_velocity).all():
            ERROR(f"No value for {str_date}")
            continue
        else:
            # plot dots at gps sites location
            fig.plot(x=df.x, y=df.y, style="c0.1c", fill="white", pen="black")
            # plot red dots at gps sites location with missing data
            fig.plot(x=df.x[~np.isnan(df.east_velocity)], y=df.y[~np.isnan(df.east_velocity)], style="c0.1c", fill="red", pen="black")
            # plot velocity vectors

        arrow_shape = ".4c+n50/0.8+gblack+e+z+qk"
        # plot dots at gps sites location
        fig.plot(x=df.x, y=df.y, style="c0.1c", fill="white", pen="black")
        # plot red dots at gps sites location with missing data
        fig.plot(x=df.x[~np.isnan(df.east_velocity)], y=df.y[~np.isnan(df.east_velocity)], style="c0.1c", fill="red", pen="black")
        # plot velocity vectors
        fig.velo(data=df, region=region, pen="1p,black", uncertaintyfill="lightblue1", line=False, spec=arrow_scale, projection=projection, vector=arrow_shape, no_clip=no_clip)
        # plot scale vector
        arrow_scale_scale = f"e{gps_scale:.3f}/0./12"
        fig.velo(data=df_scale, region=region, pen="1p,black", uncertaintyfill="lightblue1", line=False, spec=arrow_scale_scale, projection=projection, vector=arrow_shape, no_clip=no_clip)
        # plot earthquake location is provided
        if eq_lon_lat is not None:
            fig.plot(x=eq_lon_lat[0], y=eq_lon_lat[1], style="a0.5c", fill="yellow", pen="black")

        # show
        if show in ['notebook', 'external']:
            fig.show(method='external')

        # save
        if save_dir is not None:
            try:
                fig.savefig(f"{save_dir}/ivel_{str_date}.pdf")
            except Exception as e:
                ERROR(f"Could not save {save_dir}/ivel_{str_date}.pdf. Error: {e}")
