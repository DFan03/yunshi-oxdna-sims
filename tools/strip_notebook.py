#!/usr/bin/env python3
"""Strip outputs, execution counts and volatile metadata from Jupyter notebooks (in place or to --out).

Usage: python tools/strip_notebook.py nb1.ipynb [nb2.ipynb ...]
Only code/markdown cell sources are kept, so notebooks stay tiny and diff cleanly.
"""
import json
import sys


def strip(nb):
    for c in nb.get("cells", []):
        if c.get("cell_type") == "code":
            c["outputs"] = []
            c["execution_count"] = None
        c["metadata"] = {}
        c.get("attachments") and c.pop("attachments")
    ks = nb.get("metadata", {}).get("kernelspec")
    nb["metadata"] = {"kernelspec": ks} if ks else {}
    return nb


def main(paths):
    for p in paths:
        with open(p, encoding="utf-8") as f:
            nb = json.load(f)
        with open(p, "w", encoding="utf-8") as f:
            json.dump(strip(nb), f, indent=1, ensure_ascii=False)
            f.write("\n")
        print("stripped", p)


if __name__ == "__main__":
    main(sys.argv[1:])
