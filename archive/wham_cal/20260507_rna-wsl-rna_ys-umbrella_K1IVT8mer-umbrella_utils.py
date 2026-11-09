#!/usr/bin/env python3
"""
wham_com_completed.py

Run from inside the umbrella directory (where window_00, window_01, ... live).

What it does:
  - For each window_* that has a non-empty prod/trajectory.dat ("completed window"):
      * Reads prod/forces.txt  -> r0, stiff, anchor_a, anchor_b
      * Reads prod/trajectory.dat frame-by-frame:
          - anchor-anchor distance  (umbrella / WHAM coordinate)
          - hydrogen bond count     between group A and group B
                                    (the two lines of kl.csv)
  - Performs 1D WHAM over the anchor-anchor coordinate.
  - Outputs:
      * wham_profile.dat      : x_sim, x_nm, F(kBT), P
      * wham_profile.png      : PMF (left y-axis) + mean HB (right y-axis)
      * hb_distributions.png  : per-window HB count distributions (violin plot)
      * hb_per_window.dat     : window, mean dist (nm), mean/std/median/min/max HB

HB detection
------------
kl.csv must contain exactly two lines, each a comma-separated list of
nucleotide indices (0-based).  Line 1 -> group A, line 2 -> group B.

A hydrogen bond between nucleotide i (group A) and j (group B) is counted when:
    base-base distance          < HB_DIST_CUTOFF  (sim units)
    a1_i · (-a1_j)              > HB_COS_CUTOFF   (antiparallel pairing vectors)

Bonds within the same group are NOT counted.

oxDNA trajectory line format:
    cx cy cz  a1x a1y a1z  a3x a3y a3z  ...
"""

import os
import re
import glob
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# -----------------------------
# CONFIG
# -----------------------------

WINDOW_PATTERN  = "window_*"
PROD_SUBDIR     = "prod"
TRAJ_FILENAME   = "trajectory.dat"
FORCES_FILENAME = "forces.txt"

# HB calculation
KL_CSV          = "kl.csv"   # two-line CSV: line 1 = group A, line 2 = group B

# HB geometric cutoffs (oxDNA sim units / dimensionless)
HB_COS_CUTOFF   = 0.7        # cos(θ) threshold for a1_i · (-a1_j)
HB_DIST_CUTOFF  = 2.0        # max base-base distance (sim units, ~1.7 nm)

# WHAM / histogram parameters
NBINS = 80

# Sampling controls
FRAME_STRIDE          = 1
MAX_FRAMES_PER_WINDOW = 10000

# oxDNA distance unit -> nm
SIM_TO_NM = 0.8518

BETA = 1.0

# Temperature for ΔG reporting
TEMPERATURE_K = 293.0
R_KJMOLK      = 8.314e-3
RT_KJMOL      = R_KJMOLK * TEMPERATURE_K

# Output files
OUT_DATA        = "wham_profile.dat"
OUT_PNG         = "wham_profile.png"
OUT_HB_DIST_PNG = "hb_distributions.png"
OUT_HB_DAT      = "hb_per_window.dat"


# -----------------------------
# Helpers
# -----------------------------

def load_kl_two_groups(path):
    """
    Load kl.csv with exactly two non-empty lines.
    Returns (group_a, group_b) as two np.int arrays.
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
            f"kl.csv must contain exactly 2 non-empty lines of indices, "
            f"found {len(groups)}."
        )

    print(f"[load_kl_two_groups] Group A: {len(groups[0])} nucleotides")
    print(f"[load_kl_two_groups] Group B: {len(groups[1])} nucleotides")
    return groups[0], groups[1]


def get_r0_stiff_anchors(forces_path):
    """Parse forces.txt -> (r0, stiff, anchor_a, anchor_b)."""
    print(f"[get_r0_stiff_anchors] Reading {forces_path}")
    r0 = stiff = anchor_a = anchor_b = None
    with open(forces_path) as f:
        for line in f:
            s = line.strip().replace(" ", "")
            if   s.startswith("r0=")           and r0       is None:
                try: r0       = float(s.split("=")[1])
                except ValueError: pass
            elif s.startswith("stiff=")        and stiff    is None:
                try: stiff    = float(s.split("=")[1])
                except ValueError: pass
            elif s.startswith("particle=")     and anchor_a is None:
                try: anchor_a = int(s.split("=")[1])
                except ValueError: pass
            elif s.startswith("ref_particle=") and anchor_b is None:
                try: anchor_b = int(s.split("=")[1])
                except ValueError: pass
            if all(v is not None for v in [r0, stiff, anchor_a, anchor_b]):
                break
    missing = [n for n, v in [("r0", r0), ("stiff", stiff),
                               ("particle", anchor_a), ("ref_particle", anchor_b)]
               if v is None]
    if missing:
        raise ValueError(f"Could not find {missing} in {forces_path}")
    print(f"  r0={r0}  stiff={stiff}  anchor_a={anchor_a}  anchor_b={anchor_b}")
    return r0, stiff, anchor_a, anchor_b


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
    """
    Count hydrogen bonds between group_a and group_b only (inter-group).

    For every pair (i in group_a, j in group_b):
      - base-base distance < dist_cut
      - a1_i · (-a1_j)    > cos_cut
    """
    pos_a = positions[group_a]   # (Na, 3)
    pos_b = positions[group_b]   # (Nb, 3)
    a1_a  = a1_vecs[group_a]     # (Na, 3)
    a1_b  = a1_vecs[group_b]     # (Nb, 3)

    # Pairwise distances: (Na, Nb)
    diff = pos_a[:, None, :] - pos_b[None, :, :]
    dist = np.linalg.norm(diff, axis=-1)

    # Pairwise cos(θ) = a1_i · (-a1_j): (Na, Nb)
    dot = -np.einsum('ik,jk->ij', a1_a, a1_b)

    close  = dist < dist_cut
    paired = dot  > cos_cut

    return int(np.sum(close & paired))


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


# -----------------------------
# WHAM
# -----------------------------

def wham_1d(x_by_window, r0_by_window, k_by_window, nbins):
    all_x      = np.concatenate(x_by_window)
    xmin, xmax = all_x.min(), all_x.max()
    print(f"[wham_1d] x range (sim): {xmin:.4f} .. {xmax:.4f}")

    bins        = np.linspace(xmin, xmax, nbins + 1)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    dx          = bins[1] - bins[0]
    n_windows   = len(x_by_window)

    hist = np.zeros((n_windows, nbins))
    N_i  = np.zeros(n_windows)
    for i in range(n_windows):
        h, _ = np.histogram(x_by_window[i], bins=bins)
        hist[i] = h
        N_i[i]  = x_by_window[i].size

    U = 0.5 * k_by_window[:, None] * (bin_centers[None, :] - r0_by_window[:, None])**2
    f = np.zeros(n_windows)

    print("[wham_1d] Iterating...")
    for it in range(10000):
        denom = (N_i[:, None] * np.exp(BETA * f[:, None] - BETA * U)).sum(axis=0)
        P     = np.zeros(nbins)
        mask  = denom > 0
        P[mask] = hist.sum(axis=0)[mask] / denom[mask]
        norm = (P * dx).sum()
        if norm > 0:
            P /= norm

        f_new = np.array([
            -np.log(v) if (v := (P * np.exp(-BETA * U[i]) * dx).sum()) > 0 else f[i]
            for i in range(n_windows)
        ])
        f_new -= f_new.min()
        delta  = np.max(np.abs(f_new - f))
        f      = f_new

        if (it + 1) % 100 == 0:
            print(f"  iter {it+1}: max|Δf| = {delta:.3e}")
        if delta < 1e-6:
            print(f"[wham_1d] Converged at iter {it+1}")
            break
    else:
        print("[wham_1d] WARNING: did not converge")

    F        = np.full(nbins, np.nan)
    pos_mask = P > 0
    F[pos_mask] = -np.log(P[pos_mask])
    F -= np.nanmin(F)
    return bin_centers, F, P


# -----------------------------
# Plotting
# -----------------------------

def plot_wham_and_hb(bin_centers_nm, F, mean_x_nm, mean_hb, win_labels):
    """PMF + mean HB overlay on twin axes."""
    fig, ax1 = plt.subplots(figsize=(9, 5))

    color_pmf = "steelblue"
    ax1.plot(bin_centers_nm, F, "-", color=color_pmf, linewidth=2, label="PMF (WHAM)")
    ax1.set_xlabel("Anchor–anchor distance (nm)")
    ax1.set_ylabel("Free energy ($k_B T$, min = 0)", color=color_pmf)
    ax1.tick_params(axis="y", labelcolor=color_pmf)

    ax2      = ax1.twinx()
    color_hb = "tomato"
    valid    = ~np.isnan(mean_hb)
    if valid.any():
        ax2.plot(mean_x_nm[valid], mean_hb[valid], "o--",
                 color=color_hb, linewidth=1.5, markersize=6, label="Mean HB (A↔B)")
    ax2.set_ylabel("Mean inter-group hydrogen bonds", color=color_hb)
    ax2.tick_params(axis="y", labelcolor=color_hb)

    l1, lb1 = ax1.get_legend_handles_labels()
    l2, lb2 = ax2.get_legend_handles_labels()
    ax1.legend(l1 + l2, lb1 + lb2, loc="upper right")

    plt.title("WHAM PMF and mean inter-group hydrogen bonds per window")
    fig.tight_layout()
    plt.savefig(OUT_PNG, dpi=300)
    plt.close()
    print(f"[plot] Saved {OUT_PNG}")


def plot_hb_distributions(hb_by_window, mean_x_nm, win_labels):
    """
    Violin plot of per-frame inter-group HB distributions.
    One violin per window; x-axis shows window name + mean distance (nm).
    """
    valid_data   = []
    valid_x      = []
    valid_labels = []

    for hb_arr, x_nm, lbl in zip(hb_by_window, mean_x_nm, win_labels):
        if hb_arr is not None and len(hb_arr) > 0:
            valid_data.append(hb_arr)
            valid_x.append(x_nm)
            valid_labels.append(lbl)

    if not valid_data:
        print("[plot_hb_distributions] No HB data to plot.")
        return

    n   = len(valid_data)
    fig, ax = plt.subplots(figsize=(max(8, n * 0.7), 5))

    parts = ax.violinplot(valid_data, positions=range(n),
                          showmedians=True, showextrema=True)

    for pc in parts["bodies"]:
        pc.set_facecolor("tomato")
        pc.set_alpha(0.6)
    for key in ("cmedians", "cmins", "cmaxes", "cbars"):
        if key in parts:
            parts[key].set_color("darkred")
            parts[key].set_linewidth(1.2)

    means = [arr.mean() for arr in valid_data]
    ax.scatter(range(n), means, color="black", zorder=3, s=30, label="Mean")

    tick_labels = [f"{lbl}\n{x:.2f} nm" for lbl, x in zip(valid_labels, valid_x)]
    ax.set_xticks(range(n))
    ax.set_xticklabels(tick_labels, fontsize=8)

    ax.set_ylabel("Inter-group hydrogen bond count per frame")
    ax.set_xlabel("Window  (mean anchor distance)")
    ax.set_title("Inter-group HB count distribution per umbrella window  (A↔B only)")
    ax.legend(loc="upper right")
    ax.yaxis.grid(True, alpha=0.3)

    fig.tight_layout()
    plt.savefig(OUT_HB_DIST_PNG, dpi=300)
    plt.close()
    print(f"[plot] Saved {OUT_HB_DIST_PNG}")


# -----------------------------
# Main
# -----------------------------

def main():
    print("=== wham_com_completed.py (v5) ===")
    print(f"CWD         : {os.getcwd()}")
    print(f"Temperature : {TEMPERATURE_K} K  (RT = {RT_KJMOL:.3f} kJ/mol)")

    group_a, group_b = load_kl_two_groups(KL_CSV)

    completed_windows = find_completed_windows(WINDOW_PATTERN, PROD_SUBDIR, TRAJ_FILENAME)
    if not completed_windows:
        print("ERROR: No completed windows found.")
        return

    print(f"\n[main] {len(completed_windows)} completed window(s).")

    x_by_window  = []
    r0_by_window = []
    k_by_window  = []
    win_indices  = []
    mean_x_sim   = []
    hb_by_window = []
    win_labels   = []

    for w in completed_windows:
        prod_dir    = os.path.join(w, PROD_SUBDIR)
        traj_file   = os.path.join(prod_dir, TRAJ_FILENAME)
        forces_file = os.path.join(prod_dir, FORCES_FILENAME)

        print(f"\n[main] ── {w} ──")

        r0, stiff, anchor_a, anchor_b = get_r0_stiff_anchors(forces_file)
        r0_by_window.append(r0)
        k_by_window.append(stiff)

        try:
            w_idx = int(re.search(r'\d+', w.split("_", 1)[1]).group())
        except Exception:
            w_idx = len(win_indices)
        win_indices.append(w_idx)
        win_labels.append(w)

        distances = []
        hb_counts = []

        for frame_i, (positions, a1_vecs) in enumerate(
                parse_oxdna_traj_stream(traj_file, FRAME_STRIDE)):

            if frame_i >= MAX_FRAMES_PER_WINDOW:
                print(f"  Capped at {MAX_FRAMES_PER_WINDOW} frames.")
                break

            N = positions.shape[0]
            for label, idx in [("anchor_a", anchor_a), ("anchor_b", anchor_b),
                                ("group_a max", group_a.max()), ("group_b max", group_b.max())]:
                if idx >= N:
                    raise IndexError(f"{label} index {idx} out of range (N={N})")

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

        x_by_window.append(distances)
        hb_by_window.append(hb_counts)
        md = distances.mean()
        mean_x_sim.append(md)

        print(f"  Frames        : {len(distances)}")
        print(f"  Dist range    : {distances.min():.4f}–{distances.max():.4f} sim  |  mean: {md:.4f} ({md*SIM_TO_NM:.4f} nm)")
        print(f"  HB  mean±std  : {hb_counts.mean():.2f} ± {hb_counts.std():.2f}  |  median: {np.median(hb_counts):.1f}  |  range: {hb_counts.min()}–{hb_counts.max()}")

    if not x_by_window:
        print("ERROR: No data collected.")
        return

    r0_by_window = np.array(r0_by_window)
    k_by_window  = np.array(k_by_window)
    win_indices  = np.array(win_indices, dtype=int)
    mean_x_sim   = np.array(mean_x_sim)
    mean_x_nm    = mean_x_sim * SIM_TO_NM
    mean_hb      = np.array([arr.mean() for arr in hb_by_window])

    # WHAM
    print("\n[main] Running WHAM...")
    bin_centers_sim, F, P = wham_1d(x_by_window, r0_by_window, k_by_window, NBINS)
    bin_centers_nm = bin_centers_sim * SIM_TO_NM

    # Save WHAM profile
    print(f"[main] Writing {OUT_DATA}")
    with open(OUT_DATA, "w") as fh:
        fh.write("# x_sim  x_nm  F(kBT)  P\n")
        for xs, xn, Fi, Pi in zip(bin_centers_sim, bin_centers_nm, F, P):
            fh.write(f"{xs:.6f}  {xn:.6f}  {Fi:.6f}  {Pi:.6e}\n")

    # Save HB summary
    print(f"[main] Writing {OUT_HB_DAT}")
    with open(OUT_HB_DAT, "w") as fh:
        fh.write("# window  mean_dist_nm  mean_hb  std_hb  median_hb  min_hb  max_hb\n")
        for lbl, x_nm, hb_arr in zip(win_labels, mean_x_nm, hb_by_window):
            fh.write(f"{lbl}  {x_nm:.4f}  {hb_arr.mean():.3f}  "
                     f"{hb_arr.std():.3f}  {np.median(hb_arr):.1f}  "
                     f"{hb_arr.min()}  {hb_arr.max()}\n")

    # Plots
    plot_wham_and_hb(bin_centers_nm, F, mean_x_nm, mean_hb, win_labels)
    plot_hb_distributions(hb_by_window, mean_x_nm, win_labels)

    # ΔF
    i_first = np.argmin(win_indices)
    i_last  = np.argmax(win_indices)
    F_first = np.interp(mean_x_sim[i_first], bin_centers_sim, F)
    F_last  = np.interp(mean_x_sim[i_last],  bin_centers_sim, F)
    dF_kBT  = F_last - F_first

    print("\n[main] ==== ΔF: first → last window ====")
    print(f"  First (idx {win_indices[i_first]}): {mean_x_sim[i_first]:.4f} sim = {mean_x_nm[i_first]:.4f} nm  |  F = {F_first:.3f} k_B T  |  HB = {mean_hb[i_first]:.2f}")
    print(f"  Last  (idx {win_indices[i_last]}):  {mean_x_sim[i_last]:.4f} sim = {mean_x_nm[i_last]:.4f} nm  |  F = {F_last:.3f} k_B T  |  HB = {mean_hb[i_last]:.2f}")
    print(f"  ΔF = {dF_kBT:.3f} k_B T  ≈  {dF_kBT * RT_KJMOL:.3f} kJ/mol  (T = {TEMPERATURE_K} K)")
    print("=========================================")

    print(f"\n[main] Done.")
    print(f"  {OUT_DATA}  |  {OUT_PNG}  |  {OUT_HB_DIST_PNG}  |  {OUT_HB_DAT}")


if __name__ == "__main__":
    main()

