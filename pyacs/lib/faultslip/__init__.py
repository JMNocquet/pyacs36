"""
Various routines relating fault parameters, slip and motion across faults
"""

###############################################################################
# CUSTOM IMPORT
###############################################################################

__custom_import__ = ['v_to_rake',
                     'v_to_n_ss',
                     'strike_dip_rake_to_dir',
                     'rake_from_euler',
                     'rake_from_slip_az',
                     'slip_rake_2_ds_ss',
                     'ss_ns_2_ve_vn',
                     'geo_to_strike',
                     'unit_normal',
                     'unit_slip',
                     'fault_to_corners']

for mod in __custom_import__:
    exec('from .' + mod + ' import *')
