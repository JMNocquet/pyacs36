#!/usr/bin/env python3
"""
Generate autosummary stub .rst files for a package so that Sphinx builds
API docs with: (1) package/module description from __init__.py / module docstring,
(2) autosummary table with links to separate HTML files, (3) each linked file
contains the full function/class description. Works for both pyacs.lib.X when
X is a module or a package.
"""
from __future__ import annotations

import argparse
import inspect
import pkgutil
import sys
from pathlib import Path

# Ensure project root is on path when run from docs/tools or repo root
_script_dir = Path(__file__).resolve().parent
_repo_root = _script_dir.parent.parent
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))


def direct_children(package_name: str) -> list[str]:
    """Return full names of direct submodules/subpackages (no private)."""
    try:
        import importlib
        pkg = importlib.import_module(package_name)
    except Exception:
        # Import may fail (e.g. missing optional deps); fall back to filesystem
        return _direct_children_via_find_spec(package_name)
    if not hasattr(pkg, "__path__"):
        return []
    base_dots = package_name.count(".")
    out = []
    for m in pkgutil.iter_modules(pkg.__path__, pkg.__name__ + "."):
        if m.name.count(".") != base_dots + 1:
            continue
        if m.name.split(".")[-1].startswith("_"):
            continue
        out.append(m.name)
    if out:
        return sorted(out)
    # iter_modules sometimes returns [] (e.g. namespace packages); fall back to filesystem
    return _direct_children_via_find_spec(package_name)


def _direct_children_via_find_spec(package_name: str) -> list[str]:
    """Infer direct submodule names from package directory without importing the package."""
    import os
    import importlib.util
    try:
        spec = importlib.util.find_spec(package_name)
    except Exception:
        return []
    if not spec or not getattr(spec, "submodule_search_locations", None):
        return []
    out = []
    for search_dir in spec.submodule_search_locations:
        try:
            for name in os.listdir(search_dir):
                if name.startswith("_") or name == "__init__.py":
                    continue
                if name.endswith(".py"):
                    out.append(f"{package_name}.{name[:-3]}")
                elif os.path.isdir(os.path.join(search_dir, name)):
                    init = os.path.join(search_dir, name, "__init__.py")
                    if os.path.isfile(init):
                        out.append(f"{package_name}.{name}")
        except OSError:
            continue
    return sorted(out)


def get_public_members(fullname: str) -> list[tuple[str, str]]:
    """Return [(name, 'function'|'class'), ...] for public callables and classes."""
    try:
        import importlib
        mod = importlib.import_module(fullname)
    except Exception:
        return []
    names = []
    if hasattr(mod, "__all__"):
        names = list(mod.__all__)
    else:
        for name in dir(mod):
            if name.startswith("_"):
                continue
            obj = getattr(mod, name)
            if (
                inspect.isfunction(obj)
                or inspect.isclass(obj)
                or isinstance(obj, (staticmethod, classmethod))
            ):
                names.append(name)
    names = sorted(set(names))
    result = []
    for name in names:
        try:
            obj = getattr(mod, name)
            if inspect.isfunction(obj):
                result.append((name, "function"))
            elif inspect.isclass(obj):
                result.append((name, "class"))
            elif isinstance(obj, staticmethod):
                result.append((name, "function"))
            elif isinstance(obj, classmethod):
                result.append((name, "function"))
            else:
                result.append((name, "function"))
        except Exception:
            result.append((name, "function"))
    return result


def stub_module_page(fullname: str, items: list[str], _toctree_prefix: str) -> str:
    """RST for a module or package page. Use :toctree: . so links resolve to api/<fullname>.html (same dir as this page)."""
    title = fullname
    underline = "=" * len(title)
    block = f"""\
{title}
{underline}

.. automodule:: {fullname}
   :no-members:
   :show-inheritance:

"""
    if items:
        # :toctree: . = same directory as this doc → docnames api/pyacs.lib.* so HTML is api/pyacs.lib.*.html
        # ~ prefix makes Sphinx display short name (e.g. cal2datetime) instead of full path (pyacs.lib.astrotime.cal2datetime)
        items_block = "\n".join(f"   ~{x}" for x in items)
        block += f"""\
.. autosummary::
   :toctree: .
   :nosignatures:

{items_block}
"""
    return block


def stub_leaf_module_page(fullname: str) -> str:
    """RST for a leaf module: one page with full member docs so api/pyacs.lib.X.html has function description."""
    title = fullname
    underline = "=" * len(title)
    return f"""\
{title}
{underline}

.. automodule:: {fullname}
   :members:
   :undoc-members:
   :show-inheritance:
"""


def stub_member_page(parent_module: str, member_name: str, kind: str) -> str:
    """RST for a single function/class page (separate HTML file with full description)."""
    title = f"{parent_module}.{member_name}"
    underline = "=" * len(title)
    directive = "autofunction" if kind == "function" else "autoclass"
    return f"""\
{title}
{underline}

.. currentmodule:: {parent_module}

.. {directive}:: {member_name}
"""


def generate_stubs(
    package_name: str,
    out_dir: Path,
    toctree_prefix: str,
    written: set[str] | None = None,
) -> None:
    """Generate module/package pages (description + autosummary) and separate member pages (full description)."""
    if written is None:
        written = set()
    children = direct_children(package_name)

    if not children:
        # Leaf module: single page with full member docs (automodule :members:) so api/pyacs.lib.X.html has function description
        rst_path = out_dir / f"{package_name}.rst"
        if package_name not in written:
            rst_path.parent.mkdir(parents=True, exist_ok=True)
            rst_path.write_text(
                stub_leaf_module_page(package_name),
                encoding="utf-8",
            )
            written.add(package_name)
        return

    # Package: description + autosummary of submodules (full names); :toctree: . → api/pyacs.lib.*.html
    rst_path = out_dir / f"{package_name}.rst"
    if package_name not in written:
        rst_path.parent.mkdir(parents=True, exist_ok=True)
        rst_path.write_text(
            stub_module_page(package_name, children, toctree_prefix),
            encoding="utf-8",
        )
        written.add(package_name)
    for child in children:
        generate_stubs(child, out_dir, toctree_prefix, written)


def main() -> None:
    ap = argparse.ArgumentParser(description="Generate autosummary stub RST files for a package")
    ap.add_argument("--pkg", default="pyacs.lib", help="Package name (e.g. pyacs.lib)")
    ap.add_argument("--out", required=True, help="Output directory (e.g. docs/source/api/_generated/lib)")
    ap.add_argument(
        "--toctree-prefix",
        default="api/_generated/lib",
        help="Docname prefix for toctree entries (e.g. api/_generated/lib)",
    )
    args = ap.parse_args()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    generate_stubs(args.pkg, out_dir, args.toctree_prefix)
    print(f"Generated stubs for {args.pkg} under {out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
