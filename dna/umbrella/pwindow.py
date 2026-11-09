#!/usr/bin/env python3
"""
plot_window.py

Usage:
    python plot_window.py WIN_ID

Examples:
    python plot_window.py 00
    python plot_window.py 3

Run from inside the umbrella directory (where window_00, window_01, ... live).

What it does for the chosen window:
  - Reads:
      window_XX/prod/trajectory.dat
      window_XX/prod/forces.txt
      kl.csv  (two lines: line 1 = group A, line 2 = group B)
  - For each frame:
      * anchor-anchor distance  (biased coordinate)
      * inter-group HB count    (A↔B only)
  - Plots both vs frame index on a dual-axis figure:
      left  y-axis : anchor-anchor distance (nm) + r0 target line
      right y-axis : inter-group HB count
  - Saves PNG: window_XX_distances.png
"""

import os
import sys
import re
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# -----------------------------
# CONFIG
# -----------------------------
KL_CSV    = "kl.csv"    # two-line CSV in the umbrella pwd
SIM_TO_NM = 0.8518

# HB geometric cutoffs
HB_COS_CUTOFF  = 0.7
HB_DIST_CUTOFF = 2.0


# -----------------------------
# Helpers
# -----------------------------

def load_kl_two_groups(path):
    """Load kl.csv -> (group_a, group_b) as np.int arrays."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"kl.csv not found: {path}")
    groups = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            tokens = re.split(r"[,\s]+", line)
            indices = [int(t) for t in tokens if t.strip().lstrip("-").isdigit()]
            if indices:
                groups.append(np.array(indices, dtype=int))
    if len(groups) != 2:
        raise ValueError(f"kl.csv must have exactly 2 non-empty lines, found {len(groups)}.")
    print(f"[load_kl_two_groups] Group A: {len(groups[0])}  Group B: {len(groups[1])} nucleotides")
    return groups[0], groups[1]


def get_anchors_and_r0(forces_path):
    """Parse forces.txt -> (anchor_a, anchor_b, r0)."""
    if not os.path.isfile(forces_path):
        raise FileNotFoundError(f"forces.txt not found: {forces_path}")
    anchor_a = anchor_b = r0 = None
    with open(forces_path) as f:
        for line in f:
            s = line.strip().replace(" ", "")
            if s.startswith("particle=") and anchor_a is None:
                try:    anchor_a = int(s.split("=")[1])
                except ValueError: pass
            elif s.startswith("ref_particle=") and anchor_b is None:
                try:    anchor_b = int(s.split("=")[1])
                except ValueError: pass
            elif s.startswith("r0=") and r0 is None:
                try:    r0 = float(s.split("=")[1])
                except ValueError: pass
            if all(v is not None for v in [anchor_a, anchor_b, r0]):
                break
    missing = [n for n, v in [("particle", anchor_a), ("ref_particle", anchor_b), ("r0", r0)]
               if v is None]
    if missing:
        raise ValueError(f"Could not find {missing} in {forces_path}")
    print(f"  anchor_a={anchor_a}  anchor_b={anchor_b}  r0={r0:.4f} sim ({r0*SIM_TO_NM:.4f} nm)")
    return anchor_a, anchor_b, r0


def parse_oxdna_traj_stream(filepath):
    """
    Yield (positions, a1_vectors) per frame as np.array([N,3]) each.
    oxDNA line: cx cy cz  a1x a1y a1z  a3x a3y a3z  ...
    """
    with open(filepath) as f:
        in_frame  = False
        box_next  = False
        pos_rows  = []
        a1_rows   = []

        for line in f:
            if line.startswith("t ="):
                if in_frame and pos_rows:
                    yield (np.array(pos_rows, dtype=float),
                           np.array(a1_rows,  dtype=float))
                in_frame = True
                box_next = True
                pos_rows = []
                a1_rows  = []
                continue
            if not in_frame:
                continue
            if box_next:
                box_next = False
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                pos_rows.append([float(parts[0]), float(parts[1]), float(parts[2])])
                a1_rows.append( [float(parts[3]), float(parts[4]), float(parts[5])])
            except ValueError:
                continue

        if in_frame and pos_rows:
            yield (np.array(pos_rows, dtype=float),
                   np.array(a1_rows,  dtype=float))


def count_hb_frame(positions, a1_vecs, group_a, group_b, cos_cut, dist_cut):
    """Count inter-group HBs (i ∈ group_a, j ∈ group_b) for one frame."""
    pos_a = positions[group_a]
    pos_b = positions[group_b]
    a1_a  = a1_vecs[group_a]
    a1_b  = a1_vecs[group_b]
    diff  = pos_a[:, None, :] - pos_b[None, :, :]
    dist  = np.linalg.norm(diff, axis=-1)
    dot   = -np.einsum('ik,jk->ij', a1_a, a1_b)
    return int(np.sum((dist < dist_cut) & (dot > cos_cut)))


# -----------------------------
# Main
# -----------------------------

def main():
    if len(sys.argv) != 2:
        print("Usage: python plot_window.py WIN_ID")
        print("Example: python plot_window.py 00  or  python plot_window.py 3")
        sys.exit(1)

    win_str = sys.argv[1]
    try:
        win_index = int(win_str)
    except ValueError:
        print(f"ERROR: WIN_ID must be an integer, got {win_str!r}")
        sys.exit(1)

    win_name    = f"window_{win_index:02d}"
    prod_dir    = os.path.join(os.getcwd(), win_name, "prod")
    traj_file   = os.path.join(prod_dir, "trajectory.dat")
    forces_file = os.path.join(prod_dir, "forces.txt")

    print(f"=== plot_window.py  {win_name} ===")
    print(f"CWD: {os.getcwd()}")

    for label, path in [("prod dir",    prod_dir),
                        ("trajectory",  traj_file),
                        ("forces.txt",  forces_file)]:
        if not os.path.exists(path):
            print(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    group_a, group_b = load_kl_two_groups(KL_CSV)
    anchor_a, anchor_b, r0_sim = get_anchors_and_r0(forces_file)
    r0_nm = r0_sim * SIM_TO_NM

    anchor_dist_nm = []
    hb_counts      = []

    for frame_i, (positions, a1_vecs) in enumerate(parse_oxdna_traj_stream(traj_file)):
        N = positions.shape[0]
        for label, idx in [("anchor_a", anchor_a), ("anchor_b", anchor_b),
                            ("group_a max", group_a.max()), ("group_b max", group_b.max())]:
            if idx >= N:
                raise IndexError(f"{label}={idx} out of range (N={N}) in frame {frame_i}")

        anchor_dist_nm.append(
            float(np.linalg.norm(positions[anchor_a] - positions[anchor_b])) * SIM_TO_NM
        )
        hb_counts.append(
            count_hb_frame(positions, a1_vecs, group_a, group_b,
                           HB_COS_CUTOFF, HB_DIST_CUTOFF)
        )

        if (frame_i + 1) % 100 == 0:
            print(f"  {frame_i+1} frames processed...")

    if not anchor_dist_nm:
        print("ERROR: No frames found in trajectory.")
        sys.exit(1)

    anchor_dist_nm = np.array(anchor_dist_nm)
    hb_counts      = np.array(hb_counts, dtype=int)
    frames         = np.arange(len(anchor_dist_nm))

    print(f"[main] Frames: {len(frames)}")
    print(f"  Dist (nm) : {anchor_dist_nm.min():.4f} – {anchor_dist_nm.max():.4f}  mean={anchor_dist_nm.mean():.4f}")
    print(f"  HB mean±std: {hb_counts.mean():.2f} ± {hb_counts.std():.2f}  range={hb_counts.min()}–{hb_counts.max()}")

    # ── Plot ──────────────────────────────────────────────────────────────
    fig, ax1 = plt.subplots(figsize=(11, 5))

    color_dist = "steelblue"
    color_hb   = "tomato"

    # Left axis: anchor distance
    ax1.plot(frames, anchor_dist_nm, color=color_dist, alpha=0.85,
             linewidth=0.8, label="Anchor–anchor distance")
    ax1.axhline(r0_nm, color="navy", linestyle="--", linewidth=1.2,
                label=f"r0 target = {r0_nm:.2f} nm")
    ax1.set_xlabel("Frame index")
    ax1.set_ylabel("Anchor–anchor distance (nm)", color=color_dist)
    ax1.tick_params(axis="y", labelcolor=color_dist)

    # Right axis: HB count
    ax2 = ax1.twinx()
    ax2.plot(frames, hb_counts, color=color_hb, alpha=0.7,
             linewidth=0.8, label="Inter-group HB (A↔B)")
    ax2.set_ylabel("Inter-group HB count", color=color_hb)
    ax2.tick_params(axis="y", labelcolor=color_hb)
    ax2.yaxis.set_major_locator(plt.MaxNLocator(integer=True))

    # Combined legend
    l1, lb1 = ax1.get_legend_handles_labels()
    l2, lb2 = ax2.get_legend_handles_labels()
    ax1.legend(l1 + l2, lb1 + lb2, loc="upper right", fontsize=8)

    ax1.set_title(f"{win_name}: anchor distance and inter-group HB vs frame")
    ax1.grid(True, alpha=0.25)

    fig.tight_layout()
    out_png = os.path.join(os.getcwd(), f"{win_name}_distances.png")
    print(f"[main] Saving {out_png}")
    plt.savefig(out_png, dpi=300)
    plt.close()
    print("[main] Done.")


if __name__ == "__main__":
    main()

