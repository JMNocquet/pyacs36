# L1-trend Analysis Package

This package contains modular functions for L1-trend filtering, refinement, and analysis of GPS time series data.

## Package Structure

```
pyacs/gts/lib/l1trend/
├── __init__.py              # Package initialization and imports
├── breakpoints.py           # Breakpoint detection and extraction
├── cleaning.py              # Cleaning functions for L1-trend data
├── simplification.py        # Simplification algorithms
├── optimal_refinement.py    # Optimal piecewise linear refinement
├── check_trend.py           # Quality checking and analysis
├── refinement.py            # Main refinement function
├── statistics.py            # Statistics and model evaluation
├── preprocessing.py         # Pre-processing functions
├── optimization.py          # Optimization algorithms
├── l1trendi.py             # Main L1-trend filtering function
├── optimal_l1trend_workflow.py  # Optimal workflow function
├── outlier_flagging.py     # Outlier detection using L1-trend
├── sum_l1_trend.py         # L1-trend summary statistics
└── README.md               # This file
```

## Functions

### `l1trend_to_breakpoints(self, tol='auto', threshold=[1.,1., 5.])`
Convert a Gts resulting from L1-trend-filtering to a dictionary of breakpoints.

### `clean_l1trend(self, raw_gts, threshold='auto')`
Cleans breakpoints from L1-trend filtered time series.

### `simplify_l1trend(self, tolerance=.5, components='ENU')`
Remove unnecessary breakpoints from an L1-trend filtered time series.

### `optimal_pwlf_refinement(t, y, y0, yn, nbp)`
Find the optimal piecewise linear function with nbp breakpoints fitting y.

### `optimal_pwlf_refinement_fast(t, y, y0, yn)`
Fast version using l1trendi with criterion AICc.

### `check_l1_trend(ts, l1ts, component='EN', ...)`
Inspect the result from l1trend model of a time series.

### `refine_l1trend(self, rawts, lcomponent='ENU', ...)`
Main refinement function that orchestrates the entire process.

### `l1trend_optimal_workflow(self)`
Optimal workflow that combines L1-trend filtering, simplification, and refinement for best results.

### `flag_outliers_using_l1trend(self, l1trend, threshold=5)`
Flag outliers in the time series based on deviation from L1-trend representation using MAD method.

### `sum_l1_trend(self)`
Summarize L1-trend filtered time series with comprehensive statistics including breakpoints, velocities, and segment durations.

## Usage

```python
from pyacs.gts.lib.l1trend import refine_l1trend, simplify_l1trend, l1trendi

# Use the functions as methods on Gts objects
refined_ts = gts.refine_l1trend(raw_ts)
simplified_ts = gts.simplify_l1trend()
filtered_ts = gts.l1trendi(alpha='auto', criterion='BIC')

# Use the optimal workflow for best results
optimal_ts = gts.l1trend_optimal_workflow()

# Flag outliers using L1-trend
gts.flag_outliers_using_l1trend(l1trend_ts, threshold=5)

# Summarize L1-trend results
l1trend_ts.sum_l1_trend()
```

## Benefits of Modular Structure

1. **Maintainability**: Each function is in its own file, making it easier to maintain and debug
2. **Reusability**: Functions can be imported individually as needed
3. **Testing**: Each module can be tested independently
4. **Documentation**: Each file has focused documentation
5. **Performance**: Easier to optimize individual functions
6. **Collaboration**: Multiple developers can work on different modules simultaneously

## Migration from Original File

## Functions

### Core Functions
- `l1trendi()` - Main L1-trend filtering function with automatic parameter optimization
- `refine_l1trend()` - Refine L1-trend results using piecewise linear functions
- `simplify_l1trend()` - Simplify L1-trend by removing unnecessary breakpoints
- `l1trend_optimal_workflow()` - Complete optimal workflow combining all steps
- `flag_outliers_using_l1trend()` - Flag outliers based on L1-trend deviation
- `sum_l1_trend()` - Summarize L1-trend statistics and analysis

### Utility Functions
- `l1trend_to_breakpoints()` - Extract breakpoints from L1-trend results
- `clean_l1trend()` - Clean L1-trend data
- `check_l1_trend()` - Quality checking and analysis

### Optimization Functions
- `best_l1trend_golden()` - Golden section search optimization
- `best_l1trend_custom()` - Custom optimization algorithm
- `optimal_pwlf_refinement()` - Optimal piecewise linear refinement
- `optimal_pwlf_refinement_fast()` - Fast version of optimal refinement

### Preprocessing Functions
- `pre_process_test()` - Test for large offsets in time series
- `pre_process_ts()` - Pre-process time series for L1-trend

### Statistics Functions
- `get_stats_l1_model()` - Calculate statistics for model evaluation

## Refactoring

The original `refine_l1trend.py` and `l1trendi.py` files have been refactored into this modular structure. All functionality has been preserved while improving organization and maintainability.
