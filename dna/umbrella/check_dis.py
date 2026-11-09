"""
check_anchor_dist.py
--------------------
Compute the distance between two anchor particles in any oxDNA .dat file.

Usage:
    python check_anchor_dist.py <dat_file> <anchor_a> <anchor_b>

Example:
    python check_anchor_dist.py last_conf.dat 9280 9951
"""

import sys
import math


def parse_dat_coords(filepath):
    coords = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line[0] in ('t', 'b', 'E'):
                continue
            vals = line.split()
            if len(vals) < 3:
                continue
            try:
                coords.append([float(vals[0]), float(vals[1]), float(vals[2])])
            except ValueError:
                continue
    return coords


def dist(a, b):
    return math.sqrt(sum((a[k] - b[k]) ** 2 for k in range(3)))


def main():
    if len(sys.argv) != 4:
        print("Usage: python check_anchor_dist.py <dat_file> <anchor_a> <anchor_b>")
        sys.exit(1)

    dat_file = sys.argv[1]
    anchor_a = int(sys.argv[2])
    anchor_b = int(sys.argv[3])

    coords = parse_dat_coords(dat_file)

    if not coords:
        print(f"ERROR: could not parse coordinates from {dat_file}")
        sys.exit(1)

    n = len(coords)
    for idx, label in [(anchor_a, 'anchor_a'), (anchor_b, 'anchor_b')]:
        if idx >= n:
            print(f"ERROR: {label} index {idx} out of range (file has {n} particles)")
            sys.exit(1)

    pa = coords[anchor_a]
    pb = coords[anchor_b]
    d  = dist(pa, pb)

    print(f"File      : {dat_file}")
    print(f"Anchor A  : particle {anchor_a}  pos = ({pa[0]:.4f}, {pa[1]:.4f}, {pa[2]:.4f})")
    print(f"Anchor B  : particle {anchor_b}  pos = ({pb[0]:.4f}, {pb[1]:.4f}, {pb[2]:.4f})")
    print(f"Distance  : {d:.4f}  (oxDNA units)")
    print(f"          : {d * 0.8518:.4f}  nm")


if __name__ == '__main__':
    main()

