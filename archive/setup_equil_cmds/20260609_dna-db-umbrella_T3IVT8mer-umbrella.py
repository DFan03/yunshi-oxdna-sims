#!/usr/bin/env python3

from pathlib import Path

# Find all directories matching window_XX and sort them numerically
windows = sorted(
    [d for d in Path(".").iterdir() if d.is_dir() and d.name.startswith("window_")],
    key=lambda x: int(x.name.split("_")[1])
)

with open("equil_cmds", "w") as f:
    for w in windows:
        f.write(
            f"cd {w.name}/equil && "
            "/home/donglifan/oxdna/build/bin/oxDNA input.txt > run.log 2>&1 "
            "&& cd ..\n"
        )
