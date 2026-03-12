#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import sys
import re


AUTOSUMMARY_BLOCK = """\
.. autosummary::
   :toctree: {toctree}
   :nosignatures:

{items}
"""


SECTION_TITLES = ("Subpackages", "Submodules", "Subpackages and submodules")


def replace_toctree_with_autosummary(lines: list[str], toctree_dir: str) -> list[str]:
    # Find any section title we know, then the first ".. toctree::" after it.
    for title in SECTION_TITLES:
        try:
            i_sec = next(i for i, s in enumerate(lines) if s.strip() == title)
        except StopIteration:
            continue

        try:
            i_toc = next(i for i in range(i_sec, len(lines)) if lines[i].lstrip().startswith(".. toctree::"))
        except StopIteration:
            continue

        # Collect entries until next non-indented line that looks like a new section title
        entries: list[str] = []
        i = i_toc + 1
        while i < len(lines):
            s = lines[i]
            if s and (not s.startswith(" ")) and (not s.startswith("\t")):
                break
            stripped = s.strip()
            if stripped.startswith("pyacs."):
                entries.append(stripped)
            i += 1

        if not entries:
            return lines

        items = "\n".join(f"   {e}" for e in entries)
        new_block = AUTOSUMMARY_BLOCK.format(toctree=toctree_dir, items=items).splitlines()
        return lines[:i_toc] + new_block + lines[i:]

    return lines


def strip_automodule_members(lines: list[str]) -> list[str]:
    out = []
    in_automodule = False
    for s in lines:
        if s.lstrip().startswith(".. automodule::"):
            in_automodule = True
            out.append(s)
            continue

        if in_automodule:
            # automodule option lines are indented and start with ':'
            if re.match(r"^\s*:\w", s):
                if s.strip() in (":members:", ":undoc-members:"):
                    continue
                # optionally also remove :show-inheritance: if you want
                out.append(s)
                continue
            # End of automodule options when we hit a non-option (or blank then text)
            if s.strip() == "":
                out.append(s)
                continue
            in_automodule = False

        out.append(s)
    return out


def move_module_contents_to_top(lines: list[str]) -> list[str]:
    # Works if there is a "Module contents" section; if not, do nothing.
    try:
        i = next(i for i, s in enumerate(lines) if s.strip() == "Module contents")
    except StopIteration:
        return lines

    # Find end of this section: next top-level section title or end-of-file.
    start = i
    j = i + 1
    while j < len(lines):
        if lines[j] and (not lines[j].startswith(" ")) and (not lines[j].startswith("\t")):
            # likely a new section title line (followed by underline)
            if j + 1 < len(lines) and set(lines[j + 1].strip()) in ({"-"}, {"="}):
                break
        j += 1
    end = j

    module_section = lines[start:end]
    remaining = lines[:start] + lines[end:]

    # Insert after the document title (first 2 lines: title + underline)
    insert_at = 2
    if len(remaining) > 2 and remaining[2].strip() == "":
        insert_at = 3

    return remaining[:insert_at] + [""] + module_section + [""] + remaining[insert_at:]


def patch_file(rst_path: Path, toctree_dir: str) -> None:
    lines = rst_path.read_text(encoding="utf-8").splitlines()

    # 1) If there is a toctree under Submodules/Subpackages, replace by autosummary
    lines = replace_toctree_with_autosummary(lines, toctree_dir)

    # 2) Strip :members: and :undoc-members: from automodule blocks
    lines = strip_automodule_members(lines)

    # 3) Move "Module contents" section up (if present)
    lines = move_module_contents_to_top(lines)

    rst_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Patched {rst_path.name}")


def main() -> None:
    if len(sys.argv) < 2:
        print("Usage: patch_apidoc_universal.py <path/to/generated.rst> [toctree_dir]")
        sys.exit(1)
    rst = Path(sys.argv[1])
    toctree_dir = sys.argv[2] if len(sys.argv) >= 3 else "api/_generated"
    patch_file(rst, toctree_dir)


if __name__ == "__main__":
    main()