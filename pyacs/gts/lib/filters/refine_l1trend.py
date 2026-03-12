"""
L1-trend refinement functions have been moved to pyacs.gts.lib.l1trend package.

This file has been refactored into a modular structure for better organization and maintainability.

The functions have been distributed as follows:
- refine_l1trend -> pyacs.gts.lib.l1trend.refinement
- check_l1_trend -> pyacs.gts.lib.l1trend.check_trend
- l1trend_to_breakpoints -> pyacs.gts.lib.l1trend.breakpoints
- clean_l1trend -> pyacs.gts.lib.l1trend.cleaning
- simplify_l1trend -> pyacs.gts.lib.l1trend.simplification
- optimal_pwlf_refinement, optimal_pwlf_refinement_fast -> pyacs.gts.lib.l1trend.optimal_refinement

Please import from the new location:
    from pyacs.gts.lib.l1trend import refine_l1trend, check_l1_trend, l1trend_to_breakpoints, etc.
"""

# Import from new location for backward compatibility
from pyacs.gts.lib.l1trend import (
    refine_l1trend,
    check_l1_trend,
    l1trend_to_breakpoints,
    clean_l1trend,
    simplify_l1trend,
    optimal_pwlf_refinement,
    optimal_pwlf_refinement_fast
)

__all__ = [
    'refine_l1trend',
    'check_l1_trend',
    'l1trend_to_breakpoints',
    'clean_l1trend',
    'simplify_l1trend',
    'optimal_pwlf_refinement',
    'optimal_pwlf_refinement_fast'
]
