"""Super class for a collection of geodetic time series (Sgts).

Sgts holds multiple Gts instances and provides methods to operate on them
collectively.
"""

import sys
#from pyacs.gts.Gts import Gts

class Sgts:


    ###################################################################
    def __init__ (self, ts_dir='.', add_key='', verbose=True, name_filter='', read=True,sites=[],lexclude=[],type=None,xyz=True):
    ###################################################################
        """Initialize Sgts: optionally load time series from a directory.

        Parameters
        ----------
        ts_dir : str, optional
            Directory containing time series files. Default is '.'.
        add_key : str, optional
            Suffix to add to site codes.
        verbose : bool, optional
            If True, print progress. Default is True.
        name_filter : str, optional
            Filter on file/site names.
        read : bool, optional
            If True, read time series from ts_dir. Default is True.
        sites : list, optional
            If non-empty, only these sites are loaded.
        lexclude : list, optional
            Site codes to exclude.
        type : str, optional
            File type filter.
        xyz : bool, optional
            If True, load XYZ data. Default is True.
        """
        self.dir=ts_dir
        self.filter=name_filter
        self.verbose=verbose
        if read:
            self.read_ts(ts_dir=self.dir, name_filter=self.filter, add_key=add_key, verbose=verbose,sites=sites,lexclude=lexclude,type=type,xyz=xyz)



###################################################################
##  METHODS IMPORT
###################################################################

import pyacs.gts.Sgts_methods.add_offsets_dates
import pyacs.gts.Sgts_methods.append
import pyacs.gts.Sgts_methods.copy
import pyacs.gts.Sgts_methods.delts
import pyacs.gts.Sgts_methods.frame
import pyacs.gts.Sgts_methods.gts
import pyacs.gts.Sgts_methods.has_ts
import pyacs.gts.Sgts_methods.lGts
import pyacs.gts.Sgts_methods.n
import pyacs.gts.Sgts_methods.lcode
import pyacs.gts.Sgts_methods.medvel
import pyacs.gts.Sgts_methods.read_gmt
import pyacs.gts.Sgts_methods.read_gts_conf
import pyacs.gts.Sgts_methods.read_soln
import pyacs.gts.Sgts_methods.read_ts
import pyacs.gts.Sgts_methods.same_site
import pyacs.gts.Sgts_methods.save_velocity
import pyacs.gts.Sgts_methods.sel_rectangle
import pyacs.gts.Sgts_methods.sel_radius
import pyacs.gts.Sgts_methods.sel_period
import pyacs.gts.Sgts_methods.show_map
import pyacs.gts.Sgts_methods.stat_site
import pyacs.gts.Sgts_methods.sub
import pyacs.gts.Sgts_methods.to_displacement
import pyacs.gts.Sgts_methods.write_pck
import pyacs.gts.Sgts_methods.common_mode
import pyacs.gts.Sgts_methods.to_obs_tensor
import pyacs.gts.Sgts_methods.apply_coseismic
import pyacs.gts.Sgts_methods.get_unr
import pyacs.gts.Sgts_methods.info
import pyacs.gts.Sgts_methods.to_kml
import pyacs.gts.Sgts_methods.plot_data_sum
import pyacs.gts.Sgts_methods.plot_component
import pyacs.gts.Sgts_methods.dates
import pyacs.gts.Sgts_methods.sel_from_grid
import pyacs.gts.Sgts_methods.compute_common_mode_l1trend
import pyacs.gts.Sgts_methods.gts_mp
import pyacs.gts.Sgts_methods.delnone
import pyacs.gts.Sgts_methods.to_tspck
import pyacs.gts.Sgts_methods.to_tsnpz
import pyacs.gts.Sgts_methods.make_distance_matrix_from_sgts
import pyacs.gts.Sgts_methods.nearest
import pyacs.gts.Sgts_methods.show_ivel_map_gmt
import pyacs.gts.Sgts_methods.correct_offsets_from_file
import pyacs.gts.Sgts_methods.remove_observations
import pyacs.gts.Sgts_methods.get_dates
import pyacs.gts.Sgts_methods.sel_radius_eq

Sgts.add_offsets_dates  = pyacs.gts.Sgts_methods.add_offsets_dates.add_offsets_dates
Sgts.append             = pyacs.gts.Sgts_methods.append.append
Sgts.copy               = pyacs.gts.Sgts_methods.copy.copy
Sgts.delts              = pyacs.gts.Sgts_methods.delts.delts
Sgts.frame              = pyacs.gts.Sgts_methods.frame.frame
Sgts.gts                = pyacs.gts.Sgts_methods.gts.gts
Sgts.has_ts             = pyacs.gts.Sgts_methods.has_ts.has_ts
Sgts.lGts               = pyacs.gts.Sgts_methods.lGts.lGts
Sgts.n                  = pyacs.gts.Sgts_methods.n.n
Sgts.lcode              = pyacs.gts.Sgts_methods.lcode.lcode
Sgts.medvel             = pyacs.gts.Sgts_methods.medvel.medvel
Sgts.read_gmt           = pyacs.gts.Sgts_methods.read_gmt.read_gmt
Sgts.read_gts_conf      = pyacs.gts.Sgts_methods.read_gts_conf.read_gts_conf
Sgts.read_soln          = pyacs.gts.Sgts_methods.read_soln.read_soln
Sgts.read_ts            = pyacs.gts.Sgts_methods.read_ts.read_ts
Sgts.same_site          = pyacs.gts.Sgts_methods.same_site.same_site
Sgts.save_velocity      = pyacs.gts.Sgts_methods.save_velocity.save_velocity
Sgts.sel_radius         = pyacs.gts.Sgts_methods.sel_radius.sel_radius
Sgts.sel_rectangle      = pyacs.gts.Sgts_methods.sel_rectangle.sel_rectangle
Sgts.sel_period         = pyacs.gts.Sgts_methods.sel_period.sel_period
Sgts.show_map           = pyacs.gts.Sgts_methods.show_map.show_map
Sgts.stat_site          = pyacs.gts.Sgts_methods.stat_site.stat_site
Sgts.sub                = pyacs.gts.Sgts_methods.sub.sub
Sgts.to_displacement    = pyacs.gts.Sgts_methods.to_displacement.to_displacement
Sgts.write_pck          = pyacs.gts.Sgts_methods.write_pck.write_pck
Sgts.common_mode        = pyacs.gts.Sgts_methods.common_mode.common_mode
Sgts.to_obs_tensor      = pyacs.gts.Sgts_methods.to_obs_tensor.to_obs_tensor
Sgts.apply_coseismic    = pyacs.gts.Sgts_methods.apply_coseismic.apply_coseismic
Sgts.get_unr            = pyacs.gts.Sgts_methods.get_unr.get_unr
Sgts.info               = pyacs.gts.Sgts_methods.info.info
Sgts.to_kml             = pyacs.gts.Sgts_methods.to_kml.to_kml
Sgts.plot_data_sum      = pyacs.gts.Sgts_methods.plot_data_sum.plot_data_sum
Sgts.plot_component     = pyacs.gts.Sgts_methods.plot_component.plot_component
Sgts.dates              = pyacs.gts.Sgts_methods.dates.dates
Sgts.sel_from_grid      = pyacs.gts.Sgts_methods.sel_from_grid.sel_from_grid
Sgts.compute_common_mode_l1trend = pyacs.gts.Sgts_methods.compute_common_mode_l1trend.compute_common_mode_l1trend
Sgts.gts_mp             = pyacs.gts.Sgts_methods.gts_mp.gts_mp
Sgts.delnone            = pyacs.gts.Sgts_methods.delnone.delnone
Sgts.to_tsnpz           = pyacs.gts.Sgts_methods.to_tsnpz.to_tsnpz
Sgts.to_tspck           = pyacs.gts.Sgts_methods.to_tspck.to_tspck
Sgts.make_distance_matrix_from_sgts = pyacs.gts.Sgts_methods.make_distance_matrix_from_sgts.make_distance_matrix_from_sgts
Sgts.nearest           = pyacs.gts.Sgts_methods.nearest.nearest
Sgts.show_ivel_map_gmt = pyacs.gts.Sgts_methods.show_ivel_map_gmt.show_ivel_map_gmt
Sgts.correct_offsets_from_file = pyacs.gts.Sgts_methods.correct_offsets_from_file.correct_offsets_from_file
Sgts.remove_observations = pyacs.gts.Sgts_methods.remove_observations.remove_observations
Sgts.get_dates = pyacs.gts.Sgts_methods.get_dates.get_dates
Sgts.sel_radius_eq = pyacs.gts.Sgts_methods.sel_radius_eq.sel_radius_eq