#!/usr/bin/env python3
"""
Generate an ASCII report of differences between two .tex files,
with 3 lines of before/after context (with line numbers).

Usage:
    python tex_diff_ascii_report.py original.tex edited.tex report.txt
"""

import sys
import difflib
from pathlib import Path

def read_file(path: Path):
    return path.read_text(encoding="utf-8").splitlines()

def block_str(lines, i1, i2, with_nums=False):
    if with_nums:
        return "\n".join(f"{i+1:>5}: {lines[i]}" for i in range(i1, i2)).strip()
    return "\n".join(lines[i1:i2]).strip()

def context_slice(lines, idx, before=3, after=3, with_nums=False):
    start = max(0, idx - before)
    end = min(len(lines), idx + after)
    before_txt = "\n".join(
        f"{i+1:>5}: {lines[i]}" if with_nums else lines[i]
        for i in range(start, idx)
    )
    after_txt = "\n".join(
        f"{i+1:>5}: {lines[i]}" if with_nums else lines[i]
        for i in range(idx, end)
    )
    return before_txt, after_txt

def main():
    if len(sys.argv) != 4:
        print("Usage: python tex_diff_ascii_report.py original.tex edited.tex report.txt")
        sys.exit(1)

    orig_path = Path(sys.argv[1])
    edit_path = Path(sys.argv[2])
    out_path  = Path(sys.argv[3])

    orig = read_file(orig_path)
    edit = read_file(edit_path)

    sm = difflib.SequenceMatcher(None, orig, edit)
    report_lines = []
    changes = 0

    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        if tag == "equal":
            continue

        orig_block = block_str(orig, i1, i2, with_nums=True)
        edit_block = block_str(edit, j1, j2, with_nums=True)
        changes += 1

        if tag == "replace":
            header = "====== CHANGE DETECTED ======"
        elif tag == "delete":
            header = "====== DELETION DETECTED ======"
        else:  # insert
            header = "====== ADDITION DETECTED ======"

        # context
        if tag == "insert":
            before_txt = "\n".join(
                f"{i+1:>5}: {edit[i]}" for i in range(max(0, j1-3), j1)
            )
            after_txt = "\n".join(
                f"{i+1:>5}: {orig[i]}" for i in range(i1, min(len(orig), i1+3))
            )
        else:
            before_txt, after_txt = context_slice(edit, j1, before=3, after=3, with_nums=True)

        section = [
            header,
            ""
        ]

        if tag in ("replace", "delete") and orig_block:
            section += [
                "Original TEX:",
                orig_block,
                ""
            ]

        if tag in ("replace", "insert") and edit_block:
            section += [
                "Edited/Added TEX:",
                edit_block,
                ""
            ]

        section += [
            "Placement hint (raw TEX with line numbers):",
            "Before:",
            before_txt,
            "After:",
            after_txt,
            "-" * 40,
            ""
        ]

        report_lines.extend(section)

    if changes == 0:
        out_path.write_text("No differences detected between the TEX files.\n", encoding="utf-8")
        print(f"ℹ️ No differences detected. Wrote a simple report to {out_path}")
    else:
        out_path.write_text("\n".join(report_lines), encoding="utf-8")
        print(f"✅ ASCII report written to {out_path} with {changes} change(s).")

if __name__ == "__main__":
    main()
