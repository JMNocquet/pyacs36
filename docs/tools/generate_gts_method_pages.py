#!/usr/bin/env python3
"""Generate one .rst page per Gts method so :meth: links point to separate HTML files."""
from pathlib import Path

GTS_METHODS = [
    "read", "read_pos", "read_pride", "read_cats_file", "read_series",
    "write_pos", "write_cats", "to_pandas_df", "to_pytrf", "get_unr", "force_daily",
    "copy", "extract_periods", "extract_dates", "exclude_periods", "substract_ts",
    "add_obs", "xyz2neu", "neu2xyz", "reorder", "decimate", "displacement",
    "rotate", "insert_gts_data", "interpolate", "set_zero_at_date", "split", "get_coseismic",
    "detrend", "detrend_annual", "detrend_seasonal", "detrend_median", "trajectory",
    "frame", "make_model", "remove_pole", "detrend_pytrf",
    "find_offsets", "suspect_offsets", "apply_offsets", "test_offset_significance",
    "find_offsets_t_scan", "find_offsets_ivel",
    "remove_outliers", "find_outliers_simple", "find_outliers_vondrak",
    "find_outliers_sliding_window", "find_outlier_around_date",
    "realistic_sigma", "wrms", "sigma_vel_tsfit",
    "vondrak", "wiener", "spline", "median_filter", "edge_filter", "disp2vel", "ivel", "smooth",
    "l1trendi", "refine_l1trend", "clean_l1trend", "simplify_l1trend",
    "l1trend_to_breakpoints", "l1trend_optimal_workflow",
    "info", "save_velocity", "save_offsets", "plot",
]

TEMPLATE = """Gts.{method}
{underline}

.. currentmodule:: pyacs.gts.Gts

.. autoclass:: Gts
   :members: {method}
   :no-index:
"""

def main():
    script_dir = Path(__file__).resolve().parent
    source_dir = script_dir.parent / "source"
    out_dir = source_dir / "gts_methods"
    out_dir.mkdir(parents=True, exist_ok=True)
    for method in GTS_METHODS:
        rst = TEMPLATE.format(method=method, underline="=" * (len(method) + 4))
        (out_dir / f"{method}.rst").write_text(rst, encoding="utf-8")
    print(f"Generated {len(GTS_METHODS)} method pages in {out_dir}", file=__import__("sys").stderr)

if __name__ == "__main__":
    main()
