#!/usr/bin/env python3

import sys
import re

if len(sys.argv) != 3:
    print("Usage: ntremover.py forces.txt removals.txt")
    sys.exit(1)

forces_file = sys.argv[1]
removals_file = sys.argv[2]

# Load removal IDs
with open(removals_file) as f:
    removals = set(map(int, f.read().split()))

with open(forces_file) as f:
    lines = f.readlines()

output = []
block = []
inside = False

for line in lines:
    if line.strip().startswith("{"):
        inside = True
        block = [line]
        continue

    if inside:
        block.append(line)

        if line.strip().startswith("}"):
            remove = False

            for bline in block:
                m1 = re.search(r'particle\s*=\s*(\d+)', bline)
                m2 = re.search(r'ref_particle\s*=\s*(\d+)', bline)

                if m1 and int(m1.group(1)) in removals:
                    remove = True
                if m2 and int(m2.group(1)) in removals:
                    remove = True

            if not remove:
                output.extend(block)

            inside = False
            block = []

    else:
        output.append(line)

# Overwrite original file
with open(forces_file, "w") as f:
    f.writelines(output)
