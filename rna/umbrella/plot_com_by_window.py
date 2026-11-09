#!/usr/bin/env python3
"""
plot_com_by_window.py

Run from inside the umbrella directory (where window_00, window_01, ... live).

What it does:
  - Reads GROUP A / B index files (a_coord.txt, b_coord.txt).
  - Scans window_* directories.
  - For each window that has a non-empty prod/trajectory.dat:
      * Computes COM–COM distances for groups A and B over frames
      * Converts to nm
  - Plots COM distance distributions per window as horizontal boxplots:

        X-axis: distance (nm)
        Y-axis: window index (00, 01, ...)

  - Saves figure as: com_by_window.png
"""

import os
import glob
import numpy as np

# Non-interactive backend for PNG
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# -----------------------------
# CONFIG: adjust if needed
# -----------------------------
# If a_coord.txt and b_coord.txt are in the parent of umbrella/, use "../...".
# If you copied them into umbrella/, change these to "a_coord.txt" / "b_coord.txt".
GROUP_A_FILE = "../a_coord.txt"
GROUP_B_FILE = "../b_coord.txt"

WINDOW_PATTERN   = "window_*"
PROD_SUBDIR      = "prod"
TRAJ_FILENAME    = "trajectory.dat"

# Performance controls
FRAME_STRIDE          = 1     # use 1 to see full behavior; increase for speed if needed
MAX_FRAMES_PER_WINDOW = 10000 # hard cap to avoid pathological files

OUTPUT_PNG = "com_by_window.png"

SIM_TO_NM = 0.8518  # oxDNA distance unit -> nm


# -----------------------------
# Helpers
# -----------------------------

def load_indices(path):
    """Load nucleotide indices from a whitespace/CSV-like file."""
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
    Stream parser for oxDNA trajectory files.

    Yields np.array([N,3]) positions for every 'frame_stride'-th frame.
    """
    print(f"[parse_oxdna_traj_stream] Opening {filepath}")
    with open(filepath) as f:
        in_frame = False
        box_line_next = False
        frame_positions = []
        frame_idx = 0

        for line_num, line in enumerate(f, 1):
            if line.startswith("t ="):
                # Close previous frame
                if in_frame:
                    if frame_positions and (frame_idx % frame_stride == 0):
                        yield np.array(frame_positions, dtype=float)
                    frame_idx += 1

                in_frame = True
                box_line_next = True
                frame_positions = []
                continue

            if not in_frame:
                continue

            if box_line_next:
                # skip box line
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

        # Last frame
        if in_frame and frame_positions and (frame_idx % frame_stride == 0):
            yield np.array(frame_positions, dtype=float)


def com_distance(positions, idx_A, idx_B):
    """COM–COM distance for two groups of indices."""
    com_A = positions[idx_A].mean(axis=0)
    com_B = positions[idx_B].mean(axis=0)
    return float(np.linalg.norm(com_A - com_B))


def find_completed_windows(window_pattern, prod_subdir, traj_filename):
    """Return sorted list of window dirs that have a non-empty production trajectory."""
    all_windows = sorted(d for d in glob.glob(window_pattern) if os.path.isdir(d))
    completed = []

    print(f"[find_completed_windows] Scanning windows matching '{window_pattern}'")
    for w in all_windows:
        traj_path = os.path.join(w, prod_subdir, traj_filename)
        if os.path.isfile(traj_path):
            size = os.path.getsize(traj_path)
            if size > 0:
                print(f"  [OK]   {w}: {traj_path} ({size} bytes)")
                completed.append(w)
            else:
                print(f"  [SKIP] {w}: trajectory exists but is empty")
        else:
            print(f"  [SKIP] {w}: no trajectory found")

    return completed


# -----------------------------
# Main
# -----------------------------

def main():
    print("=== plot_com_by_window.py ===")
    print(f"Current directory: {os.getcwd()}")

    # Load groups
    idx_A = load_indices(GROUP_A_FILE)
    idx_B = load_indices(GROUP_B_FILE)

    # Completed windows only
    completed_windows = find_completed_windows(WINDOW_PATTERN, PROD_SUBDIR, TRAJ_FILENAME)
    if not completed_windows:
        print("ERROR: No completed windows (non-empty prod/trajectory.dat) found.")
        return

    print(f"\n[main] Using {len(completed_windows)} completed window(s):")
    for w in completed_windows:
        print(f"  - {w}")

    all_com = []
    win_ids = []

    # For each completed window, compute COM distances
    for w in completed_windows:
        prod_dir = os.path.join(w, PROD_SUBDIR)
        traj_file = os.path.join(prod_dir, TRAJ_FILENAME)

        print(f"\n[main] Processing {w}")
        print(f"  Trajectory: {traj_file}")

        distances = []
        for frame_i, positions in enumerate(parse_oxdna_traj_stream(traj_file, frame_stride=FRAME_STRIDE)):
            if frame_i >= MAX_FRAMES_PER_WINDOW:
                print(f"  Reached MAX_FRAMES_PER_WINDOW ({MAX_FRAMES_PER_WINDOW}); stopping for this window.")
                break

            n_particles = positions.shape[0]
            if idx_A.max() >= n_particles or idx_B.max() >= n_particles:
                raise IndexError(
                    f"Group indices out of range for frame {frame_i} in {traj_file} "
                    f"(N={n_particles}, max idx_A={idx_A.max()}, max idx_B={idx_B.max()})"
                )

            d = com_distance(positions, idx_A, idx_B)
            distances.append(d)

            if (frame_i + 1) % 100 == 0:
                print(f"  Processed {frame_i + 1} frames for {w}...")

        distances = np.array(distances, dtype=float)
        if distances.size == 0:
            print(f"  WARNING: No frames parsed for {w}, skipping this window.")
            continue

        # convert to nm
        distances_nm = distances * SIM_TO_NM
        all_com.append(distances_nm)

        # extract window index from name "window_XX"
        try:
            win_idx = int(w.split("_")[1])
        except (IndexError, ValueError):
            # fallback: use enumeration if parsing fails
            win_idx = len(win_ids)
        win_ids.append(win_idx)

        print(f"  Collected {len(distances_nm)} COM distances for {w}")
        print(f"  Range (nm): min={distances_nm.min():.3f}, max={distances_nm.max():.3f}")

    if not all_com:
        print("ERROR: No COM distances collected.")
        return

    # Sort windows by window index
    win_ids = np.array(win_ids, dtype=int)
    order = np.argsort(win_ids)
    win_ids_sorted = win_ids[order]
    com_sorted = [all_com[i] for i in order]

    # Plot: horizontal boxplots
    print("\n[main] Creating plot...")
    plt.figure(figsize=(10, 6))

    plt.boxplot(
        com_sorted,
        positions=win_ids_sorted,
        vert=False,          # horizontal boxes
        showfliers=False,
    )

    plt.xlabel("COM–COM distance (nm)")
    plt.ylabel("Window index")
    plt.yticks(win_ids_sorted, [f"{i:02d}" for i in win_ids_sorted])
    plt.title("COM distance ranges per completed window")
    plt.grid(axis="x", alpha=0.3)

    plt.tight_layout()
    print(f"[main] Saving figure to {OUTPUT_PNG}")
    plt.savefig(OUTPUT_PNG, dpi=300)
    plt.close()
    print("[main] Done.")


if __name__ == "__main__":
    main()
