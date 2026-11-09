#!/usr/bin/env python3
"""
plot_completed_windows.py

Run from inside the umbrella directory (where window_00, window_01, ... live).

This script:
  - Reads group A / B indices from text files.
  - Scans window_* directories.
  - For each window:
        * if prod/trajectory.dat exists and is non-empty -> treat as COMPLETED
        * else -> skip with a message
  - For each completed window:
        * reads prod/forces.txt to get r0
        * computes COM–COM distances (with stride and max frame cap)
  - Saves a PNG: umbrella_partial_windows.png
    showing only the completed windows.

Use this to tune stiffness / overlap on a subset of windows
before running everything.
"""

import os
import sys
import glob
import numpy as np

# Non-interactive backend for PNG
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# =========================
# CONFIG (edit as needed)
# =========================

# Paths to coordinate index files for the two tips.
# If a_coord.txt / b_coord.txt sit one level up: "../a_coord.txt" etc.
GROUP_A_FILE = "../a_coord.txt"
GROUP_B_FILE = "../b_coord.txt"

# Window structure (relative to where you run the script)
WINDOW_PATTERN   = "window_*"
PROD_SUBDIR      = "prod"
TRAJ_FILENAME    = "trajectory.dat"
FORCES_FILENAME  = "forces.txt"

# Performance controls
FRAME_STRIDE          = 100   # only every 100th frame → faster
MAX_FRAMES_PER_WINDOW = 500   # cap per window

# Output figure
OUTPUT_PNG = "umbrella_partial_windows.png"

# oxDNA distance unit to nm (approx for DNA model)
SIM_TO_NM = 0.8518


# =========================
# Helper functions
# =========================

def load_indices(path):
    """
    Load nucleotide indices from a whitespace/CSV-like file.
    Ignores any non-integer tokens (e.g. stray '%').
    """
    print(f"[load_indices] Loading indices from {path}")
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Index file not found: {path}")

    with open(path) as f:
        text = f.read()

    for sep in [",", ";", "\t"]:
        text = text.replace(sep, " ")

    tokens = [t for t in text.split() if t.strip()]

    indices = []
    for t in tokens:
        try:
            indices.append(int(t))
        except ValueError:
            continue

    if not indices:
        raise ValueError(f"No valid integer indices found in {path}")

    arr = np.array(indices, dtype=int)
    print(f"[load_indices] Loaded {len(arr)} indices from {path}")
    return arr


def parse_oxdna_traj_stream(filepath, frame_stride=1):
    """
    Stream parser for oxDNA trajectory/conf files.

    Yields np.array([N,3]) of positions for every 'frame_stride'-th frame.
    """
    print(f"[parse_oxdna_traj_stream] Opening {filepath}")
    with open(filepath) as f:
        in_frame = False
        box_line_next = False
        frame_positions = []
        frame_idx = 0

        for line_num, line in enumerate(f, 1):
            if line.startswith("t ="):
                # Close previous frame, if any
                if in_frame:
                    if frame_positions and (frame_idx % frame_stride == 0):
                        yield np.array(frame_positions)
                    frame_idx += 1

                in_frame = True
                box_line_next = True
                frame_positions = []
                continue

            if not in_frame:
                continue

            if box_line_next:
                # Box line; skip it
                box_line_next = False
                continue

            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
            except ValueError:
                continue
            frame_positions.append([x, y, z])

        # Emit last frame if file ended mid-frame
        if in_frame and frame_positions and (frame_idx % frame_stride == 0):
            yield np.array(frame_positions)


def com_distance(positions, idx_A, idx_B):
    """COM–COM distance for two groups of indices."""
    com_A = positions[idx_A].mean(axis=0)
    com_B = positions[idx_B].mean(axis=0)
    return float(np.linalg.norm(com_A - com_B))


def get_r0_from_forces(force_file):
    """Parse forces.txt and extract the first 'r0 = ...' value."""
    print(f"[get_r0_from_forces] Reading {force_file}")
    if not os.path.isfile(force_file):
        raise FileNotFoundError(f"Force file not found: {force_file}")
    with open(force_file) as f:
        for line in f:
            s = line.strip()
            if s.startswith("r0"):
                s = s.replace(" ", "")
                parts = s.split("=")
                if len(parts) == 2:
                    try:
                        val = float(parts[1])
                        print(f"[get_r0_from_forces] Found r0 = {val:.4f}")
                        return val
                    except ValueError:
                        pass
    raise ValueError(f"Could not find a valid 'r0 = ...' line in {force_file}")


def find_completed_windows(window_pattern, prod_subdir, traj_filename):
    """
    Find window_* dirs where prod/trajectory.dat exists and is non-empty.

    Returns list of such window directories.
    """
    all_windows = sorted(d for d in glob.glob(window_pattern) if os.path.isdir(d))
    completed = []

    print(f"[find_completed_windows] Scanning window directories matching '{window_pattern}'")
    for w in all_windows:
        traj_path = os.path.join(w, prod_subdir, traj_filename)
        if os.path.isfile(traj_path):
            size = os.path.getsize(traj_path)
            if size > 0:
                print(f"  [OK] {w}: found trajectory {traj_path} ({size} bytes)")
                completed.append(w)
            else:
                print(f"  [SKIP] {w}: trajectory exists but is empty ({traj_path})")
        else:
            print(f"  [SKIP] {w}: no trajectory at {traj_path}")

    return completed


def collect_window_data_completed(
    group_A_file,
    group_B_file,
    window_pattern,
    prod_subdir,
    traj_filename,
    forces_filename,
    frame_stride,
    max_frames_per_window,
):
    """
    Collect r0 (nm) and COM–COM distances (nm) for all COMPLETED windows
    (those with a non-empty prod/trajectory.dat).

    Returns:
      r0_nm: np.array of shape [n_completed_windows]
      distances_per_window_nm: list of 1D np.arrays
      completed_windows: list of window dir names
    """
    print("\n[collect_window_data_completed] Starting collection over completed windows...")

    idx_A = load_indices(group_A_file)
    idx_B = load_indices(group_B_file)

    completed_windows = find_completed_windows(window_pattern, prod_subdir, traj_filename)
    if not completed_windows:
        raise RuntimeError("No completed windows found (no non-empty production trajectories).")

    print(f"\n[collect_window_data_completed] Found {len(completed_windows)} completed window(s):")
    for w in completed_windows:
        print(f"  - {w}")

    r0_values = []
    all_distances = []

    for w_idx, win_dir in enumerate(completed_windows):
        print(f"\n[collect_window_data_completed] Processing completed window {w_idx}: {win_dir}")
        prod_dir = os.path.join(win_dir, prod_subdir)
        traj_file = os.path.join(prod_dir, traj_filename)
        force_file = os.path.join(prod_dir, forces_filename)

        if not os.path.isfile(force_file):
            raise FileNotFoundError(f"Missing forces.txt for completed window: {force_file}")

        r0 = get_r0_from_forces(force_file)
        r0_values.append(r0)

        print(f"[collect_window_data_completed] r0 = {r0:.4f} (sim units)")
        print(f"  Reading traj: {traj_file}")
        print(f"  frame_stride = {frame_stride}, max_frames_per_window = {max_frames_per_window}")

        distances = []
        for frame_i, frame in enumerate(parse_oxdna_traj_stream(traj_file, frame_stride=frame_stride)):
            if frame_i >= max_frames_per_window:
                print(f"  Reached max_frames_per_window ({max_frames_per_window}) → stop this window.")
                break

            n_particles = frame.shape[0]
            if idx_A.max() >= n_particles or idx_B.max() >= n_particles:
                raise IndexError(
                    f"Indices in GROUP A/B exceed number of particles ({n_particles}) "
                    f"in frame from {traj_file}"
                )

            d = com_distance(frame, idx_A, idx_B)
            distances.append(d)

            if (frame_i + 1) % 100 == 0:
                print(f"  Processed {frame_i + 1} frames for {win_dir}...")

        distances = np.array(distances, dtype=float)
        all_distances.append(distances)
        print(f"[collect_window_data_completed] Done {win_dir}, collected {len(distances)} distances.")

    r0_values = np.array(r0_values, dtype=float)
    r0_nm = r0_values * SIM_TO_NM
    distances_per_window_nm = [d * SIM_TO_NM for d in all_distances]

    print("\n[collect_window_data_completed] Finished collection for completed windows.")
    return r0_nm, distances_per_window_nm, completed_windows


def plot_umbrella_windows_partial(r0_nm, distances_per_window_nm, window_names, output_png):
    """
    Violin plot of COM–COM distance per completed window.
    """
    print("\n[plot_umbrella_windows_partial] Preparing VIOLIN plot...")

    if len(r0_nm) == 0:
        raise ValueError("No r0 values provided (empty r0_nm).")

    # Sort by r0
    order = np.argsort(r0_nm)
    r0_sorted = r0_nm[order]
    data_sorted = [distances_per_window_nm[i] for i in order]
    names_sorted = [window_names[i] for i in order]

    # Estimate width
    if len(r0_sorted) > 1:
        dr = np.min(np.diff(np.sort(r0_sorted)))
        width = dr * 0.8
    else:
        width = 0.2

    plt.figure(figsize=(10, 5))

    parts = plt.violinplot(
        data_sorted,
        positions=r0_sorted,
        widths=width,
        showmeans=True,
        showextrema=True,
        showmedians=True,
    )

    plt.xlabel("Target tip-to-tip distance r0 (nm)")
    plt.ylabel("Sampled COM–COM distance (nm)")
    plt.title("Umbrella windows (completed subset): COM–COM distance per window")
    plt.grid(axis="y", alpha=0.3)

    # Labels
    labels = [f"{name}\n{r:.1f} nm" for name, r in zip(names_sorted, r0_sorted)]
    plt.xticks(r0_sorted, labels, rotation=45, ha="right")

    plt.tight_layout()
    print(f"[plot_umbrella_windows_partial] Saving figure to {output_png}")
    plt.savefig(output_png, dpi=300)
    plt.close()
    print("[plot_umbrella_windows_partial] Plot saved.")

# =========================
# Main
# =========================

def main():
    print("=== plot_completed_windows.py ===")
    print(f"Current directory: {os.getcwd()}")
    print(f"GROUP_A_FILE: {GROUP_A_FILE}")
    print(f"GROUP_B_FILE: {GROUP_B_FILE}")
    print(f"FRAME_STRIDE = {FRAME_STRIDE}")
    print(f"MAX_FRAMES_PER_WINDOW = {MAX_FRAMES_PER_WINDOW}")

    r0_nm, distances_nm, window_names = collect_window_data_completed(
        group_A_file=GROUP_A_FILE,
        group_B_file=GROUP_B_FILE,
        window_pattern=WINDOW_PATTERN,
        prod_subdir=PROD_SUBDIR,
        traj_filename=TRAJ_FILENAME,
        forces_filename=FORCES_FILENAME,
        frame_stride=FRAME_STRIDE,
        max_frames_per_window=MAX_FRAMES_PER_WINDOW,
    )

    print("\n[main] r0 values (nm):", r0_nm)
    print("[main] Completed windows:", window_names)
    print("[main] Now plotting...")
    plot_umbrella_windows_partial(r0_nm, distances_nm, window_names, OUTPUT_PNG)
    print("\n[main] Done. Output PNG:", OUTPUT_PNG)


if __name__ == "__main__":
    main()
