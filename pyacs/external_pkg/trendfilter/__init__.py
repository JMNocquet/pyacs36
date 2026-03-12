"""Bundled trendfilter (from trend_filter-master). Used by pyacs Gts.l1trend / edge_filter."""
from trendfilter.trendfilter import trend_filter

__version__ = "0.2.1"

try:
    from trendfilter.plot_model import plot_model
except ImportError:
    plot_model = None  # bokeh not installed
try:
    from trendfilter.get_example_data import get_example_data, get_example_data_seasonal
except ImportError:
    get_example_data = None
    get_example_data_seasonal = None
