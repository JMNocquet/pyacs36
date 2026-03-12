"""
L1-trend filtering functions have been moved to pyacs.gts.lib.l1trend package.

This file has been refactored into a modular structure for better organization and maintainability.

The functions have been distributed as follows:
- get_stats_l1_model -> pyacs.gts.lib.l1trend.statistics
- pre_process_test, pre_process_ts -> pyacs.gts.lib.l1trend.preprocessing  
- best_l1trend_golden, best_l1trend_custom -> pyacs.gts.lib.l1trend.optimization
- l1trendi -> pyacs.gts.lib.l1trend.l1trendi

Please import from the new location:
    from pyacs.gts.lib.l1trend import l1trendi, get_stats_l1_model, pre_process_test, etc.
"""

# Import from new location for backward compatibility
from pyacs.gts.lib.l1trend import (
    l1trendi,
    get_stats_l1_model,
    pre_process_test,
    pre_process_ts,
    best_l1trend_golden,
    best_l1trend_custom
)

__all__ = [
    'l1trendi',
    'get_stats_l1_model', 
    'pre_process_test',
    'pre_process_ts',
    'best_l1trend_golden',
    'best_l1trend_custom'
]

