# Bundled external packages

Code in this directory is from external projects and is installed by setuptools as **top-level** packages when you `pip install pyacs`, so that existing imports (e.g. `from trendfilter import trend_filter`) keep working.

- **trendfilter**: L1 trend filtering (used by `Gts.l1trend`, `edge_filter`, etc.). Sourced from [trend_filter-master](https://github.com/dave31415/trendfilter). Depends on `cvxpy` (added to pyacs `install_requires`).
