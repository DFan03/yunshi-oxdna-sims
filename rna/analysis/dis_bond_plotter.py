#!/usr/bin/env python3

import sys
import csv
import os
import matplotlib.pyplot as plt

OUTPUT_FILE = "kl_distance_vs_KL_HBond.png" 
def safe_int_from_str(s):
    """
    Convert a string to int in a robust way:
    - allow floats like '999500000.0'
    - skip if it can't be converted
    """
    try:
        # First, try direct int
        return int(s)
    except ValueError:
        try:
            # Then try float -> int (for '999500000.0' etc.)
            return int(float(s))
        except ValueError:
            return None


def read_distance_file(path):
    """
    Read distance CSV: t, distance
    Skips header and non-numeric rows.
    """
    t_vals = []
    distances = []
    with open(path, newline="") as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            t_raw = row[0].strip()
            d_raw = row[1].strip() if len(row) > 1 else None

            t = safe_int_from_str(t_raw)
            if t is None:
                # header or bad row
                continue

            try:
                distance = float(d_raw)
            except (TypeError, ValueError):
                continue

            t_vals.append(t)
            distances.append(distance)
    return t_vals, distances


def read_bond_file(path):
    """
    Read bond CSV: t(as float-like), n_bonds, (optional third column)
    Skips header and non-numeric rows.
    """
    t_vals = []
    n_bonds = []
    with open(path, newline="") as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue

            t_raw = row[0].strip()
            b_raw = row[1].strip() if len(row) > 1 else None

            t = safe_int_from_str(t_raw)
            if t is None:
                # header or bad row
                continue

            try:
                bonds = float(b_raw)
            except (TypeError, ValueError):
                continue

            t_vals.append(t)
            n_bonds.append(bonds)
    return t_vals, n_bonds


def main():
    if len(sys.argv) != 3:
        print("Usage: dis_bond_plotter.py <distance.csv> <bond.csv>")
        sys.exit(1)

    distance_file = sys.argv[1]
    bond_file = sys.argv[2]

    # Read data
    t_dist, distances = read_distance_file(distance_file)
    t_bond, bonds = read_bond_file(bond_file)

    if not t_dist or not t_bond:
        print("No valid numeric data found in one of the input files.")
        sys.exit(1)

    # Sort by t in case input is unordered
    dist_data = sorted(zip(t_dist, distances), key=lambda x: x[0])
    bond_data = sorted(zip(t_bond, bonds), key=lambda x: x[0])

    t_dist_sorted, dist_sorted = zip(*dist_data)
    t_bond_sorted, bond_sorted = zip(*bond_data)

    # Labels from file names (without .csv)
    dist_label = os.path.splitext(os.path.basename(distance_file))[0]
    bond_label = "number of bonds between kissing loops" 

    fig, ax1 = plt.subplots(figsize=(15, 5))

    # Distance line (red)
    ax1.plot(
        t_dist_sorted,
        dist_sorted,
        color="red",
        alpha=0.5,
        linewidth=2,
        label=dist_label,
    )
    ax1.set_xlabel("timestep",size=20)
    ax1.set_ylabel("distance kissing loop",size=20, color="red")
    ax1.tick_params(axis="y", labelcolor="red")

    # Bond line on second y-axis (blue)
    ax2 = ax1.twinx()
    ax2.plot(
        t_bond_sorted,
        bond_sorted,
        color="blue",
        alpha=0.5,
        linewidth=2,
        label=bond_label,
    )
    ax2.set_ylabel("Hbonds kissing loop",size=20, color="blue")
    ax2.tick_params(axis="y", labelcolor="blue")

    # Let matplotlib autoscale both y-axes for readability
    ax1.margins(y=0.05)
    ax2.margins(y=1)
    ax2.set_ylim([-0.5,7])

    # Combine legends from both axes
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="best")

    fig.tight_layout()
    plt.savefig(OUTPUT_FILE, dpi=300, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()
