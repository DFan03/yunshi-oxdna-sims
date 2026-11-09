#!/usr/bin/env python3
"""
plot_completed_windows.py

Run from inside the umbrella directory (where window_00, window_01, ... live).

This script:
  - Reads two groups of nucleotide indices from kl.csv (two lines).
  - Scans window_*/prod/ directories for non-empty trajectory.dat files.
  - For each completed window:
        * reads prod/forces.txt  -> r0, anchor_a, anchor_b
        * reads prod/trajectory.dat frame-by-frame (strided):
            - anchor-anchor distance  (the biased coordinate)
            - inter-group HB count    (group A line 1 <-> group B line 2 of kl.csv)
  - Saves two PNGs:
        * umbrella_windows_distance.png  : horizontal boxplot of anchor-anchor distance
        * umbrella_windows_hb.png        : horizontal boxplot of inter-group HB count

HB criterion (inline, no external tools):
    base-base distance  < HB_DIST_CUTOFF  (sim units)
    a1_i · (-a1_j)      > HB_COS_CUTOFF   (antiparallel pairing vectors)
Only pairs where i ∈ group A and j ∈ group B are counted.
"""

import os
import re
import glob
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# =========================
# CONFIG
# =========================

KL_CSV          = "kl.csv"       # two-line CSV: line 1 = group A, line 2 = group B

WINDOW_PATTERN  = "window_*"
PROD_SUBDIR     = "prod"
TRAJ_FILENAME   = "trajectory.dat"
FORCES_FILENAME = "forces.txt"

# HB geometric cutoffs (oxDNA sim units / dimensionless)
HB_COS_CUTOFF   = 0.7
HB_DIST_CUTOFF  = 2.0

# Performance controls
# FRAME_STRIDE=1 reads every frame (same as plot_window.py).
# Increase only if trajectories are very large and speed is a concern.
FRAME_STRIDE          = 1
MAX_FRAMES_PER_WINDOW = 10000    # match the WHAM script cap

# oxDNA distance unit -> nm
SIM_TO_NM = 0.8518

# Output figures
OUT_DIST_PNG = "umbrella_windows_distance.png"
OUT_HB_PNG   = "umbrella_windows_hb.png"


# =========================
# Helpers
# =========================

def load_kl_two_groups(path):
    """
    Load kl.csv with exactly two non-empty lines.
    Returns (group_a, group_b) as np.int arrays.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"kl.csv not found: {path}")

    groups = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            tokens = re.split(r"[,\s]+", line)
            indices = []
            for t in tokens:
                t = t.strip()
                if t:
                    try:
                        indices.append(int(t))
                    except ValueError:
                        pass
            if indices:
                groups.append(np.array(indices, dtype=int))

    if len(groups) != 2:
        raise ValueError(
            f"kl.csv must have exactly 2 non-empty lines, found {len(groups)}."
        )
    print(f"[load_kl_two_groups] Group A: {len(groups[0])} nucleotides")
    print(f"[load_kl_two_groups] Group B: {len(groups[1])} nucleotides")
    return groups[0], groups[1]


def get_r0_and_anchors(forces_path):
    """Parse forces.txt -> (r0, anchor_a, anchor_b)."""
    if not os.path.isfile(forces_path):
        raise FileNotFoundError(f"forces.txt not found: {forces_path}")
    r0 = anchor_a = anchor_b = None
    with open(forces_path) as f:
        for line in f:
            s = line.strip().replace(" ", "")
            if s.startswith("r0=") and r0 is None:
                try:    r0 = float(s.split("=")[1])
                except ValueError: pass
            elif s.startswith("particle=") and anchor_a is None:
                try:    anchor_a = int(s.split("=")[1])
                except ValueError: pass
            elif s.startswith("ref_particle=") and anchor_b is None:
                try:    anchor_b = int(s.split("=")[1])
                except ValueError: pass
            if all(v is not None for v in [r0, anchor_a, anchor_b]):
                break
    missing = [n for n, v in [("r0", r0), ("particle", anchor_a), ("ref_particle", anchor_b)]
               if v is None]
    if missing:
        raise ValueError(f"Could not find {missing} in {forces_path}")
    print(f"  r0={r0:.4f}  anchor_a={anchor_a}  anchor_b={anchor_b}")
    return r0, anchor_a, anchor_b


def parse_oxdna_traj_stream(filepath, frame_stride=1):
    """
    Yield (positions, a1_vectors) per frame as np.array([N,3]) each.
    oxDNA line: cx cy cz  a1x a1y a1z  a3x a3y a3z  ...
    """
    with open(filepath) as f:
        in_frame  = False
        box_next  = False
        pos_rows  = []
        a1_rows   = []
        frame_idx = 0

        for line in f:
            if line.startswith("t ="):
                if in_frame and pos_rows and frame_idx % frame_stride == 0:
                    yield (np.array(pos_rows, dtype=float),
                           np.array(a1_rows,  dtype=float))
                if in_frame:
                    frame_idx += 1
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

        if in_frame and pos_rows and frame_idx % frame_stride == 0:
            yield (np.array(pos_rows, dtype=float),
                   np.array(a1_rows,  dtype=float))


def count_hb_frame(positions, a1_vecs, group_a, group_b, cos_cut, dist_cut):
    """Count inter-group HBs (i ∈ group_a, j ∈ group_b) for one frame."""
    pos_a = positions[group_a]
    pos_b = positions[group_b]
    a1_a  = a1_vecs[group_a]
    a1_b  = a1_vecs[group_b]

    diff = pos_a[:, None, :] - pos_b[None, :, :]
    dist = np.linalg.norm(diff, axis=-1)
    dot  = -np.einsum('ik,jk->ij', a1_a, a1_b)

    return int(np.sum((dist < dist_cut) & (dot > cos_cut)))


def find_completed_windows(pattern, prod_subdir, traj_filename):
    all_windows = sorted(d for d in glob.glob(pattern) if os.path.isdir(d))
    completed   = []
    print(f"[find_completed_windows] Scanning '{pattern}'")
    for w in all_windows:
        traj = os.path.join(w, prod_subdir, traj_filename)
        if os.path.isfile(traj) and os.path.getsize(traj) > 0:
            print(f"  [OK]   {w}  ({os.path.getsize(traj)} bytes)")
            completed.append(w)
        else:
            print(f"  [SKIP] {w}")
    return completed


# =========================
# Data collection
# =========================

def collect_window_data(group_a, group_b, completed_windows,
                        prod_subdir, traj_filename, forces_filename,
                        frame_stride, max_frames):
    """
    Returns:
      win_ids      : np.array of int window indices
      r0_nm        : np.array of r0 values in nm
      dist_list    : list of np.arrays of anchor-anchor distances (nm)
      hb_list      : list of np.arrays of inter-group HB counts (int)
      win_names    : list of window directory names
    """
    win_ids   = []
    r0_list   = []
    dist_list = []
    hb_list   = []
    names     = []

    for w in completed_windows:
        prod_dir    = os.path.join(w, prod_subdir)
        traj_file   = os.path.join(prod_dir, traj_filename)
        forces_file = os.path.join(prod_dir, forces_filename)

        print(f"\n[collect] ── {w} ──")
        r0, anchor_a, anchor_b = get_r0_and_anchors(forces_file)

        try:
            w_idx = int(w.split("_")[1])
        except (IndexError, ValueError):
            w_idx = len(win_ids)

        distances = []
        hb_counts = []

        for frame_i, (positions, a1_vecs) in enumerate(
                parse_oxdna_traj_stream(traj_file, frame_stride)):

            if frame_i >= max_frames:
                print(f"  Capped at {max_frames} frames.")
                break

            N = positions.shape[0]
            for label, idx in [("anchor_a", anchor_a), ("anchor_b", anchor_b),
                                ("group_a max", group_a.max()), ("group_b max", group_b.max())]:
                if idx >= N:
                    raise IndexError(f"{label}={idx} out of range (N={N})")

            distances.append(float(np.linalg.norm(
                positions[anchor_a] - positions[anchor_b])))
            hb_counts.append(count_hb_frame(
                positions, a1_vecs, group_a, group_b,
                HB_COS_CUTOFF, HB_DIST_CUTOFF))

            if (frame_i + 1) % 500 == 0:
                print(f"  {frame_i+1} frames processed...")

        distances = np.array(distances)
        hb_counts = np.array(hb_counts, dtype=int)

        if distances.size == 0:
            print(f"  WARNING: no frames parsed, skipping {w}.")
            continue

        win_ids.append(w_idx)
        r0_list.append(r0 * SIM_TO_NM)
        dist_list.append(distances * SIM_TO_NM)
        hb_list.append(hb_counts)
        names.append(w)

        print(f"  Frames     : {len(distances)}")
        print(f"  Dist (nm)  : {distances.min()*SIM_TO_NM:.3f} – {distances.max()*SIM_TO_NM:.3f}  mean={distances.mean()*SIM_TO_NM:.3f}")
        print(f"  HB mean±std: {hb_counts.mean():.2f} ± {hb_counts.std():.2f}  median={np.median(hb_counts):.1f}  range={hb_counts.min()}–{hb_counts.max()}")

    return np.array(win_ids, dtype=int), np.array(r0_list), dist_list, hb_list, names


# =========================
# Plotting
# =========================

def horizontal_boxplot(data_list, win_ids, r0_nm, win_names,
                       xlabel, title, output_png):
    """
    Horizontal boxplot — one box per window.
    Y-axis: window index.  X-axis: the measured quantity.
    """
    if not data_list:
        print(f"[horizontal_boxplot] No data, skipping {output_png}.")
        return

    order        = np.argsort(win_ids)
    ids_sorted   = win_ids[order]
    data_sorted  = [data_list[i] for i in order]
    r0_sorted    = r0_nm[order]
    names_sorted = [win_names[i] for i in order]

    fig, ax = plt.subplots(figsize=(10, max(4, len(ids_sorted) * 0.45)))

    ax.boxplot(
        data_sorted,
        positions=ids_sorted,
        vert=False,
        showfliers=False,
        patch_artist=True,
        boxprops=dict(facecolor="steelblue", alpha=0.6),
        medianprops=dict(color="navy", linewidth=1.5),
        whiskerprops=dict(color="steelblue"),
        capprops=dict(color="steelblue"),
    )

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Window index")
    ax.set_yticks(ids_sorted)
    ax.set_yticklabels(
        [f"{idx:02d}  {name}  (r0={r:.2f} nm)"
         for idx, name, r in zip(ids_sorted, names_sorted, r0_sorted)],
        fontsize=8,
    )
    ax.set_title(title)
    ax.grid(axis="x", alpha=0.3)

    fig.tight_layout()
    plt.savefig(output_png, dpi=300)
    plt.close()
    print(f"[plot] Saved {output_png}")


# =========================
# Main
# =========================

def main():
    print("=== plot_completed_windows.py (v5) ===")
    print(f"CWD         : {os.getcwd()}")
    print(f"KL_CSV      : {KL_CSV}")
    print(f"FRAME_STRIDE: {FRAME_STRIDE}  |  MAX_FRAMES: {MAX_FRAMES_PER_WINDOW}")

    group_a, group_b = load_kl_two_groups(KL_CSV)

    completed = find_completed_windows(WINDOW_PATTERN, PROD_SUBDIR, TRAJ_FILENAME)
    if not completed:
        print("ERROR: No completed windows found.")
        return

    win_ids, r0_nm, dist_nm, hb_counts, win_names = collect_window_data(
        group_a, group_b, completed,
        PROD_SUBDIR, TRAJ_FILENAME, FORCES_FILENAME,
        FRAME_STRIDE, MAX_FRAMES_PER_WINDOW,
    )

    if len(win_ids) == 0:
        print("ERROR: No data collected.")
        return

    horizontal_boxplot(
        dist_nm, win_ids, r0_nm, win_names,
        xlabel     = "Anchor–anchor distance (nm)",
        title      = "Umbrella windows — anchor-anchor distance distribution",
        output_png = OUT_DIST_PNG,
    )

    horizontal_boxplot(
        hb_counts, win_ids, r0_nm, win_names,
        xlabel     = "Inter-group HB count per frame  (A↔B)",
        title      = "Umbrella windows — inter-group hydrogen bond distribution",
        output_png = OUT_HB_PNG,
    )

    print(f"\n[main] Done.")
    print(f"  {OUT_DIST_PNG}")
    print(f"  {OUT_HB_PNG}")


if __name__ == "__main__":
    main()

