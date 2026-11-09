#!/usr/bin/env python3
"""
plot_window.py

Usage:
    python plot_window.py WIN_ID

Examples:
    python plot_window.py 00
    python plot_window.py 0

Run this from inside the umbrella directory (where window_00, window_01, ... live).

What it does for the chosen window:
  - Reads:
      window_XX/prod/trajectory.dat
      window_XX/prod/forces.txt
      a_coord.txt (group A indices, relative to parent or current dir)
      b_coord.txt (group B indices)
  - Extracts anchor indices and r0 from forces.txt
  - For each frame in trajectory.dat:
      * Computes anchor–anchor distance
      * Computes COM–COM distance between groups A and B
  - Plots distance vs frame index for both:
      * anchor distance
      * COM distance
      * horizontal line at r0
  - Saves PNG: window_XX_distances.png

Adjust GROUP_A_FILE / GROUP_B_FILE below if needed.
"""

import os
import sys
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

# oxDNA length unit -> nm (approx for DNA model)
SIM_TO_NM = 0.8518


# -----------------------------
# Helper functions
# -----------------------------

def load_indices(path):
    """
    Load nucleotide indices from a whitespace/CSV-like file.
    Ignores non-integer tokens (e.g. stray '%').
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


def parse_oxdna_traj_stream(filepath):
    """
    Stream parser for oxDNA trajectory files.

    Yields np.array([N,3]) positions for each frame.
    Each frame:
      - starts with 't =' line
      - followed by one 'box' line
      - then N lines of coordinates
    """
    print(f"[parse_oxdna_traj_stream] Opening {filepath}")
    with open(filepath) as f:
        in_frame = False
        box_line_next = False
        frame_positions = []

        for line_num, line in enumerate(f, 1):
            if line.startswith("t ="):
                # Close previous frame (if any)
                if in_frame and frame_positions:
                    yield np.array(frame_positions, dtype=float)
                # Start new frame
                in_frame = True
                box_line_next = True
                frame_positions = []
                continue

            if not in_frame:
                continue

            if box_line_next:
                # This is the box line; skip it
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
        if in_frame and frame_positions:
            yield np.array(frame_positions, dtype=float)


def get_anchors_and_r0_from_forces(forces_path):
    """
    Parse forces.txt to extract:
      - anchor_a index (first 'particle =')
      - anchor_b index (first 'ref_particle =')
      - r0 value (first 'r0 = ...')

    Returns (anchor_a, anchor_b, r0)
    """
    print(f"[get_anchors_and_r0_from_forces] Reading {forces_path}")
    if not os.path.isfile(forces_path):
        raise FileNotFoundError(f"forces.txt not found: {forces_path}")

    anchor_a = None
    anchor_b = None
    r0 = None

    with open(forces_path) as f:
        for line in f:
            s = line.strip()
            # anchors
            if s.startswith("particle"):
                # particle = 1234
                parts = s.replace(" ", "").split("=")
                if len(parts) == 2 and anchor_a is None:
                    try:
                        anchor_a = int(parts[1])
                        print(f"  Found anchor A (particle) = {anchor_a}")
                    except ValueError:
                        pass
            elif s.startswith("ref_particle"):
                # ref_particle = 5678
                parts = s.replace(" ", "").split("=")
                if len(parts) == 2 and anchor_b is None:
                    try:
                        anchor_b = int(parts[1])
                        print(f"  Found anchor B (ref_particle) = {anchor_b}")
                    except ValueError:
                        pass
            elif s.startswith("r0"):
                # r0 = ...
                parts = s.replace(" ", "").split("=")
                if len(parts) == 2 and r0 is None:
                    try:
                        r0 = float(parts[1])
                        print(f"  Found r0 = {r0}")
                    except ValueError:
                        pass

            if anchor_a is not None and anchor_b is not None and r0 is not None:
                break

    if anchor_a is None or anchor_b is None:
        raise ValueError(f"Could not find anchor indices in {forces_path}")
    if r0 is None:
        raise ValueError(f"Could not find r0 in {forces_path}")

    return anchor_a, anchor_b, r0


def com_distance(positions, idx_A, idx_B):
    """
    COM–COM distance for two groups of indices.
    positions: np.array [N, 3]
    """
    com_A = positions[idx_A].mean(axis=0)
    com_B = positions[idx_B].mean(axis=0)
    return float(np.linalg.norm(com_A - com_B))


def anchor_distance(positions, anchor_a, anchor_b):
    """
    Distance between two anchor particles.
    """
    pa = positions[anchor_a]
    pb = positions[anchor_b]
    return float(np.linalg.norm(pa - pb))


# -----------------------------
# Main
# -----------------------------

def main():
    if len(sys.argv) != 2:
        print("Usage: python plot_window.py WIN_ID")
        print("Example: python plot_window.py 00  or  python plot_window.py 0")
        sys.exit(1)

    win_str = sys.argv[1]
    try:
        win_index = int(win_str)
    except ValueError:
        print(f"ERROR: WIN_ID must be an integer-like string (e.g. '0', '00', '3', '09'), got {win_str!r}")
        sys.exit(1)

    win_name = f"window_{win_index:02d}"
    print(f"=== plot_window.py for {win_name} ===")
    print(f"Current directory: {os.getcwd()}")

    win_dir = os.path.join(os.getcwd(), win_name)
    prod_dir = os.path.join(win_dir, "prod")
    traj_file = os.path.join(prod_dir, "trajectory.dat")
    forces_file = os.path.join(prod_dir, "forces.txt")

    if not os.path.isdir(win_dir):
        print(f"ERROR: {win_name} directory not found here.")
        sys.exit(1)
    if not os.path.isdir(prod_dir):
        print(f"ERROR: {prod_dir} directory not found.")
        sys.exit(1)
    if not os.path.isfile(traj_file):
        print(f"ERROR: trajectory file not found: {traj_file}")
        sys.exit(1)
    if not os.path.isfile(forces_file):
        print(f"ERROR: forces file not found: {forces_file}")
        sys.exit(1)

    # Load group indices
    idx_A = load_indices(GROUP_A_FILE)
    idx_B = load_indices(GROUP_B_FILE)

    # Anchors and r0
    anchor_a, anchor_b, r0_sim = get_anchors_and_r0_from_forces(forces_file)
    r0_nm = r0_sim * SIM_TO_NM
    print(f"[main] anchor A = {anchor_a}, anchor B = {anchor_b}, r0 = {r0_sim:.4f} (sim units) = {r0_nm:.4f} nm")

    # Arrays for distances
    anchor_distances = []
    com_distances = []
    frame_indices = []

    # Loop over frames
    for frame_i, positions in enumerate(parse_oxdna_traj_stream(traj_file)):
        n_particles = positions.shape[0]
        if anchor_a >= n_particles or anchor_b >= n_particles:
            raise IndexError(
                f"Anchor index out of range for frame {frame_i} "
                f"(anchors {anchor_a},{anchor_b}, N={n_particles})"
            )
        if idx_A.max() >= n_particles or idx_B.max() >= n_particles:
            raise IndexError(
                f"Group indices out of range for frame {frame_i} "
                f"(max idx_A={idx_A.max()}, max idx_B={idx_B.max()}, N={n_particles})"
            )

        d_anchor = anchor_distance(positions, anchor_a, anchor_b)
        d_com = com_distance(positions, idx_A, idx_B)

        anchor_distances.append(d_anchor)
        com_distances.append(d_com)
        frame_indices.append(frame_i)

        if (frame_i + 1) % 100 == 0:
            print(f"  Processed {frame_i + 1} frames...")

    anchor_distances = np.array(anchor_distances, dtype=float)
    com_distances = np.array(com_distances, dtype=float)
    frame_indices = np.array(frame_indices, dtype=int)

    if len(frame_indices) == 0:
        print("ERROR: No frames found in trajectory.")
        sys.exit(1)

    print(f"[main] Total frames processed: {len(frame_indices)}")

    # Convert to nm
    anchor_nm = anchor_distances * SIM_TO_NM
    com_nm = com_distances * SIM_TO_NM

    # Plot
    fig, ax = plt.subplots(figsize=(10, 5))

    ax.plot(frame_indices, anchor_nm, label="Anchor–anchor distance (nm)", alpha=0.8)
    ax.plot(frame_indices, com_nm, label="COM–COM distance (nm)", alpha=0.8)
    ax.axhline(r0_nm, color="k", linestyle="--", label=f"r0 (anchor target) = {r0_nm:.2f} nm")

    ax.set_xlabel("Frame index")
    ax.set_ylabel("Distance (nm)")
    ax.set_title(f"{win_name}: anchor and COM distances vs frame")
    ax.grid(True, alpha=0.3)
    ax.legend()

    plt.tight_layout()

    out_png = os.path.join(os.getcwd(), f"{win_name}_distances.png")
    print(f"[main] Saving plot to {out_png}")
    plt.savefig(out_png, dpi=300)
    plt.close()
    print("[main] Done.")


if __name__ == "__main__":
    main()
