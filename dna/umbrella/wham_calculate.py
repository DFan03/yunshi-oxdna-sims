#!/usr/bin/env python3
"""
wham_com_completed.py

Run from inside the umbrella directory (where window_00, window_01, ... live).

What it does:
  - Reads group A / B indices (COM groups) from a_coord.txt and b_coord.txt.
  - For each window_* that has a non-empty prod/trajectory.dat ("completed window"):
      * Reads prod/forces.txt to get:
          - r0 (target distance, sim units)
          - stiff (spring constant, dimensionless in oxDNA units)
      * Reads prod/trajectory.dat and computes COM–COM distance for groups A & B.
  - Treats COM distance as the umbrella coordinate with bias:
        U_i(x) = 0.5 * stiff_i * (x - r0_i)^2   (dimensionless, β = 1)
  - Performs 1D WHAM over COM coordinate using only completed windows.
  - Outputs:
      * com_wham_profile.dat : columns (x_sim, x_nm, F, P)
        where F is in k_B T units (min shifted to 0).
      * com_wham_profile.png : PMF F(x) vs x (nm).
  - Additionally:
      * Computes mean COM distance for the lowest-index window and highest-index window.
      * Interpolates F at those two mean COM distances.
      * Prints ΔF (last - first) in k_B T and kJ/mol (assuming TEMPERATURE_K).
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
FORCES_FILENAME  = "forces.txt"

# WHAM / histogram parameters
NBINS = 80

# Sampling controls
FRAME_STRIDE          = 1      # use 1 for full data; increase (>1) to speed up
MAX_FRAMES_PER_WINDOW = 10000  # hard cap per window for safety

# oxDNA distance unit -> nm
SIM_TO_NM = 0.8518

# In oxDNA, energies are in k_B T units, so β = 1
BETA = 1.0

# Temperature for reporting ΔG in kJ/mol
TEMPERATURE_K = 293.0
R_KJMOLK = 8.314e-3   # kJ/mol/K
RT_KJMOL = R_KJMOLK * TEMPERATURE_K

# Output files
OUT_DATA = "com_wham_profile.dat"
OUT_PNG  = "com_wham_profile.png"


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
                # Close previous frame, if any
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
                # Skip box line
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
            yield np.array(frame_positions, dtype=float)


def com_distance(positions, idx_A, idx_B):
    """COM–COM distance for two groups of indices."""
    com_A = positions[idx_A].mean(axis=0)
    com_B = positions[idx_B].mean(axis=0)
    return float(np.linalg.norm(com_A - com_B))


def get_r0_and_stiff_from_forces(forces_path):
    """
    Parse forces.txt to extract r0 and stiff from the first mutual_trap block.
    Returns (r0, stiff) in simulation units.
    """
    print(f"[get_r0_and_stiff_from_forces] Reading {forces_path}")
    if not os.path.isfile(forces_path):
        raise FileNotFoundError(f"forces.txt not found: {forces_path}")

    r0 = None
    stiff = None

    with open(forces_path) as f:
        for line in f:
            s = line.strip()
            if s.startswith("r0"):
                parts = s.replace(" ", "").split("=")
                if len(parts) == 2 and r0 is None:
                    try:
                        r0 = float(parts[1])
                        print(f"  Found r0 = {r0}")
                    except ValueError:
                        pass
            elif s.startswith("stiff"):
                parts = s.replace(" ", "").split("=")
                if len(parts) == 2 and stiff is None:
                    try:
                        stiff = float(parts[1])
                        print(f"  Found stiff = {stiff}")
                    except ValueError:
                        pass
            if r0 is not None and stiff is not None:
                break

    if r0 is None or stiff is None:
        raise ValueError(f"Could not find r0 and stiff in {forces_path}")

    return r0, stiff


def find_completed_windows(window_pattern, prod_subdir, traj_filename):
    """
    Return sorted list of window dirs that have a non-empty production trajectory.
    """
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
# WHAM implementation (1D)
# -----------------------------

def wham_1d(x_by_window, r0_by_window, k_by_window, nbins):
    """
    Simple 1D WHAM implementation.

    x_by_window: list of 1D arrays [x_ij] (sim units) for each window i
    r0_by_window: array [r0_i]
    k_by_window : array [k_i] (stiffness)
    nbins: number of bins for x

    Returns:
      bin_centers (sim units),
      F (free energy, dimensionless, min shifted to 0),
      P (normalized probability)
    """
    n_windows = len(x_by_window)
    assert n_windows == len(r0_by_window) == len(k_by_window)

    # Determine global x range
    all_x = np.concatenate(x_by_window)
    xmin = all_x.min()
    xmax = all_x.max()
    print(f"[wham_1d] Global x range (sim units): {xmin:.4f} .. {xmax:.4f}")

    # Create bins
    bins = np.linspace(xmin, xmax, nbins + 1)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    dx = bins[1] - bins[0]

    # Histograms and counts
    hist = np.zeros((n_windows, nbins), dtype=float)
    N_i = np.zeros(n_windows, dtype=float)

    for i in range(n_windows):
        h, _ = np.histogram(x_by_window[i], bins=bins)
        hist[i, :] = h
        N_i[i] = x_by_window[i].size
        print(f"[wham_1d] Window {i}: N_i = {N_i[i]}, non-empty bins = {(h > 0).sum()}")

    # Bias potentials U_i(x_k) at bin centers
    # U_i(x) = 0.5 * k_i * (x - r0_i)^2
    U = np.zeros((n_windows, nbins), dtype=float)
    for i in range(n_windows):
        U[i, :] = 0.5 * k_by_window[i] * (bin_centers - r0_by_window[i])**2

    # Initialize f_i
    f = np.zeros(n_windows, dtype=float)

    # WHAM iteration
    max_iter = 10000
    tol = 1e-6

    print("[wham_1d] Starting WHAM iterations...")
    for it in range(max_iter):
        # Compute P_k
        numerator = hist.sum(axis=0)  # sum_i n_i(x_k)
        # denominator: sum_i N_i * exp(β f_i - β U_i(x_k))
        denom = np.zeros(nbins, dtype=float)
        for i in range(n_windows):
            denom += N_i[i] * np.exp(BETA * f[i] - BETA * U[i, :])

        # Avoid division by zero
        mask = denom > 0
        P = np.zeros(nbins, dtype=float)
        P[mask] = numerator[mask] / denom[mask]

        # Normalize P
        norm = (P * dx).sum()
        if norm > 0:
            P /= norm

        # Update f_i
        f_new = np.zeros_like(f)
        for i in range(n_windows):
            # e^{-f_i} = ∑_k P_k e^{-U_i(x_k)} Δx
            integrand = P * np.exp(-BETA * U[i, :])
            val = (integrand * dx).sum()
            if val > 0:
                f_new[i] = -np.log(val)
            else:
                f_new[i] = f[i]

        # Shift f to avoid drift
        f_new -= f_new.min()

        delta = np.max(np.abs(f_new - f))
        f = f_new

        if (it + 1) % 100 == 0:
            print(f"  iter {it+1}: max |Δf| = {delta:.3e}")
        if delta < tol:
            print(f"[wham_1d] Converged in {it+1} iterations (max |Δf| = {delta:.3e})")
            break
    else:
        print("[wham_1d] WARNING: did not converge within max_iter")

    # Final F(x) = -log P(x) + const; set min(F)=0
    F = np.full_like(P, np.nan)
    pos = P > 0
    F[pos] = -np.log(P[pos])
    F -= np.nanmin(F)

    return bin_centers, F, P


# -----------------------------
# Main
# -----------------------------

def main():
    print("=== wham_com_completed.py ===")
    print(f"Current directory: {os.getcwd()}")
    print(f"Assumed temperature for ΔG: {TEMPERATURE_K} K (RT ≈ {RT_KJMOL:.3f} kJ/mol)")

    # Load COM groups
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

    x_by_window = []
    r0_by_window = []
    k_by_window  = []

    win_indices = []      # numeric indices (0, 1, 2, ...)
    mean_x_sim   = []     # mean COM distance per window (sim units)

    # Collect COM data and umbrella parameters
    for w in completed_windows:
        prod_dir = os.path.join(w, PROD_SUBDIR)
        traj_file = os.path.join(prod_dir, TRAJ_FILENAME)
        forces_file = os.path.join(prod_dir, FORCES_FILENAME)

        print(f"\n[main] Processing {w}")
        print(f"  Trajectory: {traj_file}")
        print(f"  Forces    : {forces_file}")

        if not os.path.isfile(forces_file):
            raise FileNotFoundError(f"Missing forces.txt for window {w}: {forces_file}")

        r0, stiff = get_r0_and_stiff_from_forces(forces_file)
        r0_by_window.append(r0)
        k_by_window.append(stiff)

        # Parse window index from name "window_XX"
        try:
            w_idx = int(w.split("_")[1])
        except (IndexError, ValueError):
            w_idx = len(win_indices)
        win_indices.append(w_idx)

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

            if (frame_i + 1) % 500 == 0:
                print(f"  Processed {frame_i + 1} frames for {w}...")

        distances = np.array(distances, dtype=float)
        if distances.size == 0:
            print(f"  WARNING: No frames parsed for {w}, skipping this window.")
            continue

        x_by_window.append(distances)
        mean_d = distances.mean()
        mean_x_sim.append(mean_d)

        print(f"  Collected {len(distances)} COM distances (sim units).")
        print(f"  Range (sim units): min={distances.min():.3f}, max={distances.max():.3f}")
        print(f"  Mean COM (sim units): {mean_d:.3f} ({mean_d*SIM_TO_NM:.3f} nm)")

    if not x_by_window:
        print("ERROR: No COM distances collected from any completed window.")
        return

    r0_by_window = np.array(r0_by_window, dtype=float)
    k_by_window  = np.array(k_by_window, dtype=float)
    win_indices  = np.array(win_indices, dtype=int)
    mean_x_sim   = np.array(mean_x_sim, dtype=float)

    # WHAM
    print("\n[main] Running 1D WHAM on COM distances (sim units)...")
    bin_centers_sim, F, P = wham_1d(x_by_window, r0_by_window, k_by_window, NBINS)

    # Convert bin centers to nm
    bin_centers_nm = bin_centers_sim * SIM_TO_NM

    # Save data
    print(f"[main] Writing WHAM profile to {OUT_DATA}")
    with open(OUT_DATA, "w") as f:
        f.write("# x_sim  x_nm  F  P\n")
        for x_s, x_n, Fi, Pi in zip(bin_centers_sim, bin_centers_nm, F, P):
            f.write(f"{x_s:.6f}  {x_n:.6f}  {Fi:.6f}  {Pi:.6e}\n")

    # Plot PMF
    print(f"[main] Plotting PMF to {OUT_PNG}")
    plt.figure(figsize=(8, 5))
    plt.plot(bin_centers_nm, F, "-o", markersize=3)
    plt.xlabel("COM–COM distance (nm)")
    plt.ylabel("Free energy F(x) (k_B T; shifted min=0)")
    plt.title("WHAM PMF along COM distance (completed windows)")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=300)
    plt.close()

    # -------------------------
    # Automatic ΔF between first and last window
    # -------------------------
    # Identify first and last window by index
    i_first = np.argmin(win_indices)
    i_last  = np.argmax(win_indices)

    x_first_sim = mean_x_sim[i_first]
    x_last_sim  = mean_x_sim[i_last]
    x_first_nm  = x_first_sim * SIM_TO_NM
    x_last_nm   = x_last_sim * SIM_TO_NM

    # Interpolate F at these positions
    F_first = np.interp(x_first_sim, bin_centers_sim, F)
    F_last  = np.interp(x_last_sim,  bin_centers_sim, F)
    dF_kBT  = F_last - F_first
    dG_kJ   = dF_kBT * RT_KJMOL

    print("\n[main] ==== Automatic ΔF between first and last window (by index) ====")
    print(f"  First window index: {win_indices[i_first]}  (mean COM ≈ {x_first_sim:.3f} sim units = {x_first_nm:.3f} nm)")
    print(f"  Last  window index: {win_indices[i_last]}  (mean COM ≈ {x_last_sim:.3f} sim units = {x_last_nm:.3f} nm)")
    print(f"  F(first) ≈ {F_first:.3f} k_B T")
    print(f"  F(last)  ≈ {F_last:.3f} k_B T")
    print(f"  ΔF = F(last) - F(first) ≈ {dF_kBT:.3f} k_B T")
    print(f"     ≈ {dG_kJ:.3f} kJ/mol  (assuming T = {TEMPERATURE_K} K)")
    print("==================================================================")


    print("[main] Done.")
    print(f"  Data  : {OUT_DATA}")
    print(f"  Figure: {OUT_PNG}")


if __name__ == "__main__":
    main()
