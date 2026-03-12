#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
import re

API_DIR = Path("docs/source/api")

AUTOMODULE_BLOCK = """
.. automodule:: {modname}
   :members:
   :undoc-members:
   :show-inheritance:
"""

def is_package_index(rst: str) -> bool:
    # apidoc package pages usually contain a "Subpackages and submodules" section and a toctree
    return ("Subpackages and submodules" in rst) and (".. toctree::" in rst)

def has_automodule_for(modname: str, rst: str) -> bool:
    return f".. automodule:: {modname}" in rst

def main() -> None:
    for p in API_DIR.glob("pyacs.*.rst"):
        modname = p.stem  # e.g. "pyacs.lib.astrotime"
        rst = p.read_text(encoding="utf-8")

        if is_package_index(rst) and not has_automodule_for(modname, rst):
            rst = rst.rstrip() + AUTOMODULE_BLOCK.format(modname=modname)
            p.write_text(rst, encoding="utf-8")

if __name__ == "__main__":
    main()
