"""
    Super Class of Geodetic Time Series Class & methods (Sgts)
    Sgts is a record of Gts and enables to apply methods at the same time to various Gts
        
"""

import sys
from pyacs.gts.Gts import Gts

class Sgts:


    ###################################################################
    def __init__ (self, ts_dir='.', add_key='', verbose=True, name_filter='', read=True,sites=[],lexclude=[],type=None,xyz=True):
    ###################################################################
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

#
# 
#     ###################################################################
#     def stack(self , verbose=True):
#     ###################################################################
#         """
#         Compute the stack of the Gts in the current Sgts
#          
#         :return Gts: the stacked Gts with '_STK' code
#         
