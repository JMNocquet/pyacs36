"""
L1-trend analysis and refinement tools for pyacs.

This package contains functions for L1-trend filtering, refinement, and analysis
of GPS time series data.
"""

from .breakpoints import l1trend_to_breakpoints
from .cleaning import clean_l1trend
from .simplification import simplify_l1trend, simplify_l1trend_with_fisher_test
from .refinement import refine_l1trend
from .optimal_refinement import optimal_pwlf_refinement, optimal_pwlf_refinement_fast, preprocess_timeseries_for_optimization
from .check_trend import check_l1_trend
from .statistics import get_stats_l1_model
from .preprocessing import pre_process_test, pre_process_ts
from .optimization import best_l1trend_golden, best_l1trend_custom
from .l1trendi import l1trendi
from .optimal_l1trend_workflow import l1trend_optimal_workflow
from .outlier_flagging import flag_outliers_using_l1trend
from .sum_l1_trend import sum_l1_trend
from .l1trend2d import l1_trendfilter_2d, find_optimal_lambda, plot_lambda_selection

__all__ = [
    'l1trend_to_breakpoints',
    'clean_l1trend', 
    'simplify_l1trend',
    'simplify_l1trend_with_fisher_test',
    'refine_l1trend',
    'optimal_pwlf_refinement',
    'optimal_pwlf_refinement_fast',
    'preprocess_timeseries_for_optimization',
    'check_l1_trend',
    'get_stats_l1_model',
    'pre_process_test',
    'pre_process_ts',
    'best_l1trend_golden',
    'best_l1trend_custom',
    'l1trendi',
    'l1trend_optimal_workflow',
    'flag_outliers_using_l1trend',
    'sum_l1_trend',
    'l1_trendfilter_2d',
    'find_optimal_lambda',
    'plot_lambda_selection'
]
