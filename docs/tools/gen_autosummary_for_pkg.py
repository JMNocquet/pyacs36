#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib
import pkgutil
import sys
from pathlib import Path

# Ensure project root is on path when run from docs/tools or repo root
_script_dir = Path(__file__).resolve().parent
_repo_root = _script_dir.parent.parent
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))


TEMPLATE = """\
{title}
{underline}

.. currentmodule:: {pkg}

.. autosummary::
   :toctree: {toctree}
   :recursive:
   :nosignatures:

{items}
"""


def list_direct_children(pkg_name: str) -> list[str]:
    pkg = importlib.import_module(pkg_name)
    if not hasattr(pkg, "__path__"):
        return []
    base_dots = pkg_name.count(".")
    out: list[str] = []
    for m in pkgutil.iter_modules(pkg.__path__, pkg.__name__ + "."):
        # direct children only: pyacs.lib.<child>
        if m.name.count(".") == base_dots + 1:
            out.append(m.name.split(".")[-1])
    # exclude private modules
    out = [x for x in out if not x.startswith("_")]
    return sorted(out)


def list_direct_children_full_names(pkg_name: str) -> list[str]:
    """Return full module names (e.g. pyacs.lib.astrotime) for direct children."""
    pkg = importlib.import_module(pkg_name)
    if not hasattr(pkg, "__path__"):
        return []
    base_dots = pkg_name.count(".")
    out: list[str] = []
    for m in pkgutil.iter_modules(pkg.__path__, pkg.__name__ + "."):
        if m.name.count(".") == base_dots + 1:
            if m.name.split(".")[-1].startswith("_"):
                continue
            out.append(m.name)
    return sorted(out)


TOCTREE_TEMPLATE = """\
{title}
{underline}

.. toctree::
   :maxdepth: 2

{entries}
"""


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pkg", required=True, help="e.g. pyacs.lib")
    ap.add_argument("--out", required=True, help="e.g. docs/source/libraries.rst")
    ap.add_argument("--title", required=True)
    ap.add_argument("--toctree", required=True, help="path relative to docs/source, e.g. api/_generated/lib")
    ap.add_argument(
        "--format",
        choices=("autosummary", "toctree"),
        default="toctree",
        help="Output autosummary (Sphinx generates stubs) or toctree (use pre-generated stubs). Default: toctree",
    )
    args = ap.parse_args()

    if args.format == "toctree":
        full_names = list_direct_children_full_names(args.pkg)
        entries = "\n".join(f"   {args.toctree.rstrip('/')}/{name}" for name in full_names)
        rst = TOCTREE_TEMPLATE.format(
            title=args.title,
            underline="=" * len(args.title),
            entries=entries or "   # (no submodules)",
        )
    else:
        items = list_direct_children(args.pkg)
        items_block = "\n".join(f"   {m}" for m in items) if items else "   # (no submodules found)"
        rst = TEMPLATE.format(
            title=args.title,
            underline="=" * len(args.title),
            pkg=args.pkg,
            toctree=args.toctree,
            items=items_block,
        )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(rst, encoding="utf-8")


if __name__ == "__main__":
    main()
