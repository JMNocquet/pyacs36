"""
Offset detection and correction for Gts time series.

This package distributes the former single-module offset functionality
into separate modules. All public Gts methods are re-exported here for
backward compatibility.

Use::

    import pyacs.gts.lib.offset
    Gts.apply_offsets = pyacs.gts.lib.offset.apply_offsets
    # or
    from pyacs.gts.lib.offset import apply_offsets, find_offsets
"""

from .apply_offsets import apply_offsets
from .delete_small_offsets import delete_small_offsets
from .estimate_local_offset import estimate_local_offset
from .find_offsets import find_offsets
from .find_offsets_ivel import find_offsets_ivel
from .find_offsets_t_scan import find_offsets_t_scan
from .find_time_offsets import find_time_offsets
from .local_offset_robust import local_offset_robust
from .suspect_offsets import suspect_offsets
from .suspect_offsets_mf import suspect_offsets_mf
from .test_offset_significance import test_offset_significance
from .test_offsets import test_offsets

__all__ = [
    'apply_offsets',
    'delete_small_offsets',
    'estimate_local_offset',
    'find_offsets',
    'find_offsets_ivel',
    'find_offsets_t_scan',
    'find_time_offsets',
    'local_offset_robust',
    'suspect_offsets',
    'suspect_offsets_mf',
    'test_offset_significance',
    'test_offsets',
]
