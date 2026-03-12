#!/usr/bin/env bash
# Build PYACS API documentation with Sphinx (autosummary only; no sphinx-apidoc).
# Each HTML page has: (1) package/module description from __init__.py / docstring,
# (2) autosummary table with links to separate HTML files that contain the full
# function/class description. Works for both pyacs.lib.X when X is a module or a package.

set -e
cd "$(dirname "$0")"

rm -rf docs/source/api/_generated
rm -rf docs/build/html/*

# Generate libraries.rst with toctree pointing at api/pyacs.lib.* (docnames under api/)
python3 docs/tools/gen_autosummary_for_pkg.py \
  --pkg pyacs.lib \
  --out docs/source/libraries.rst \
  --title "PYACS core libraries" \
  --toctree api \
  --format toctree

# Generate stub .rst files under docs/source/api/ so docnames are api/pyacs.lib.* and
# built HTML is api/pyacs.lib.*.html with full function/class descriptions
python3 docs/tools/generate_autosummary_stubs.py \
  --pkg pyacs.lib \
  --out docs/source/api \
  --toctree-prefix api

# conf.py has autosummary_generate = False so these stubs are not overwritten
python3 -m sphinx -E -b html docs/source docs/build/html
