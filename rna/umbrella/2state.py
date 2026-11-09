"""
2state.py

Compute a two-state (open / closed) free energy from:

  - an existing WHAM profile (x_sim, x_nm, F(kBT), P)
  - per-frame hydrogen bond counts from umbrella trajectories

States:
  - "open"   : HB == 0
  - "closed" : HB == 6
Anything in between (1–5 HBs) is not counted as a state in the final
two-state populations (they influence P(open|r), P(closed|r) but are
removed by renormalization).

Run from inside the umbrella directory (where window_* live).

Inputs:
  - WHAM profile file (e.g. com_wham_profile.dat), via --wham-file
  - Optional distance window [xmin, xmax] in nm, via --xmin-nm / --xmax-nm
    Distances outside this range are ignored and not included in the output.

Outputs:
  - open_closed_profile.dat :
        x_sim, x_nm, P_r, P_open_given_r, P_closed_given_r,
        N_tot, N_open, N_closed
  - open_closed_profile.png :
        PMF (F) + P(open|r) and P(closed|r) vs distance
  - Summary of P_open, P_closed, ΔF in kBT and kJ/mol (printed to stdout)
"""

import os
import re
import glob
import argparse
from collections import namedtuple

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# OAT / oxpy imports
import oxpy
from oxDNA_analysis_tools.UTILS.oat_multiprocesser import oat_multiprocesser, get_chunk_size
from oxDNA_analysis_tools.UTILS.RyeReader import describe, get_input_parameter

# -----------------------------
# CONFIG (can be edited)
# -----------------------------

WINDOW_PATTERN  = "window_*"
PROD_SUBDIR     = "prod"
TRAJ_FILENAME   = "trajectory.dat"
FORCES_FILENAME = "forces.txt"
INPUT_FILENAME  = "input.txt"
KL_CSV          = "kl.csv"

HB_NCPUS        = 1
FRAME_STRIDE          = 1
MAX_FRAMES_PER_WINDOW = 10000

SIM_TO_NM = 0.8518

TEMPERATURE_K = 293.0
R_KJMOLK      = 8.314e-3
RT_KJMOL      = R_KJMOLK * TEMPERATURE_K

OUT_OPEN_CLOSED_DAT = "open_closed_profile.dat"
OUT_OPEN_CLOSED_PNG = "open_closed_profile.png"

OPEN_HB   = 0
CLOSED_HB = 6


# -----------------------------
# Helpers: general
# -----------------------------

def load_kl_indices(path):
    """Load flat list of nucleotide indices from kl.csv (any number of lines)."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"kl.csv not found: {path}")
    with open(path) as f:
        content = f.read()
    for sep in [",", ";", "\t", "\n"]:
        content = content.replace(sep, " ")
    indices = []
    for tok in content.split():
        try:
            indices.append(int(tok))
        except ValueError:
            continue
    if not indices:
        raise ValueError(f"No valid indices in {path}")
    print(f"[load_kl_indices] {len(indices)} indices loaded from {path}")
    return indices


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
    """Yield np.array([N,3]) positions for every frame_stride-th frame."""
    with open(filepath) as f:
        in_frame  = False
        box_next  = False
        pos_rows  = []
        frame_idx = 0
        for line in f:
            if line.startswith("t ="):
                if in_frame and pos_rows and frame_idx % frame_stride == 0:
                    yield np.array(pos_rows, dtype=float)
                if in_frame:
                    frame_idx += 1
                in_frame = True
                box_next = True
                pos_rows = []
                continue
            if not in_frame:
                continue
            if box_next:
                box_next = False
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                pos_rows.append([float(parts[0]), float(parts[1]), float(parts[2])])
            except ValueError:
                continue
        if in_frame and pos_rows and frame_idx % frame_stride == 0:
            yield np.array(pos_rows, dtype=float)


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
# OAT HB calculation
# -----------------------------

ComputeContext = namedtuple(
    "ComputeContext", ["traj_info", "top_info", "input_file", "indices"]
)


def _oat_compute(ctx: ComputeContext, chunk_size: int, chunk_id: int):
    """
    Uses oxpy hb_list observable to compute per-frame HB counts.

    Only bonds where BOTH nucleotides are in ctx.indices are counted.
    """
    indices_set = set(ctx.indices)

    # prod_dir is the directory that contains input.txt
    prod_dir  = os.path.dirname(os.path.abspath(ctx.input_file))
    saved_cwd = os.getcwd()

    try:
        os.chdir(prod_dir)

        with oxpy.Context():
            inp = oxpy.InputFile()
            inp.init_from_filename(ctx.input_file)
            inp["log_file"]               = "/dev/null"
            inp["list_type"]              = "cells"
            inp["trajectory_file"]        = ctx.traj_info.path
            inp["analysis_bytes_to_skip"] = str(
                ctx.traj_info.idxs[chunk_id * chunk_size].offset
            )
            inp["confs_to_analyse"] = str(chunk_size)
            inp["analysis_data_output_1"] = (
                "{ \n"
                " name = stdout \n"
                " print_every = 1e10 \n"
                " col_1 = { \n"
                "  id = my_obs \n"
                "  type = hb_list \n"
                " } \n"
                "} "
            )

            backend = oxpy.analysis.AnalysisBackend(inp)
            frame_total_bonds = []

            while backend.read_next_configuration():
                step = backend.config_info().current_step
                raw  = (
                    backend.config_info()
                    .get_observable_by_id("my_obs")
                    .get_output_string(step)
                )
                lines       = raw.strip().split("\n")
                pairs_lines = lines[1:]
                total_bonds = 0
                for p in pairs_lines:
                    p = p.strip()
                    if not p:
                        continue
                    cols = p.split()
                    try:
                        a, b = int(cols[0]), int(cols[1])
                    except (ValueError, IndexError):
                        continue
                    if a in indices_set and b in indices_set:
                        total_bonds += 1
                frame_total_bonds.append(total_bonds)

    finally:
        os.chdir(saved_cwd)

    return np.array(frame_total_bonds, dtype=int)


def compute_hb_oat(prod_dir, input_filename, traj_filename, indices, ncpus=1):
    """
    Run OAT hb_list over a production trajectory.

    Returns np.array of per-frame total_bonds (int), or None on failure.
    """
    prod_dir   = os.path.abspath(prod_dir)
    input_file = os.path.join(prod_dir, input_filename)
    traj_file  = os.path.join(prod_dir, traj_filename)

    if not os.path.isfile(input_file):
        print(f"  [HB-OAT] WARNING: {input_file} not found, skipping HB.")
        return None

    # Resolve topology relative to prod_dir (it may be a bare filename in input.txt)
    top_raw  = get_input_parameter(input_file, "topology")
    top_file = top_raw if os.path.isabs(top_raw) else os.path.join(prod_dir, top_raw)

    if not os.path.isfile(top_file):
        print(f"  [HB-OAT] WARNING: topology not found at {top_file}, skipping HB.")
        return None

    top_info, traj_info = describe(top_file, traj_file)

    ctx        = ComputeContext(traj_info, top_info, input_file, indices)
    chunk_size = get_chunk_size()
    n_confs    = traj_info.nconfs

    total_bonds = np.empty(n_confs, dtype=int)

    def callback(i, r):
        start = i * chunk_size
        end   = start + len(r)
        total_bonds[start:end] = r

    oat_multiprocesser(n_confs, ncpus, _oat_compute, callback, ctx)

    print(f"  [HB-OAT] {n_confs} frames  |  mean={total_bonds.mean():.2f}  "
          f"std={total_bonds.std():.2f}  "
          f"range={total_bonds.min()}–{total_bonds.max()}")
    return total_bonds


# -----------------------------
# WHAM loading (no recomputation)
# -----------------------------

def load_wham_profile(path):
    """
    Load existing WHAM profile: columns x_sim, x_nm, F(kBT), P.

    Returns:
      x_sim : np.array
      x_nm  : np.array
      F     : np.array
      P     : np.array
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"WHAM profile file not found: {path}")

    xs_list, xn_list, F_list, P_list = [], [], [], []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            xs, xn, Fi, Pi = parts[:4]
            xs_list.append(float(xs))
            xn_list.append(float(xn))
            F_list.append(float(Fi))
            P_list.append(float(Pi))

    x_sim = np.array(xs_list)
    x_nm  = np.array(xn_list)
    F     = np.array(F_list)
    P     = np.array(P_list)

    print(f"[load_wham_profile] Loaded {len(x_sim)} bins from {path}")
    return x_sim, x_nm, F, P


def construct_bin_edges_from_centers(bin_centers):
    """
    Given bin centers c_j, construct bin edges e_k (len = N+1) such that:
      e_0     = c_0 - 0.5*(c_1 - c_0)
      e_k     = 0.5*(c_{k-1} + c_k) for 1 <= k <= N-1
      e_N     = c_{N-1} + 0.5*(c_{N-1} - c_{N-2})
    """
    c = np.asarray(bin_centers)
    if c.size < 2:
        raise ValueError("Need at least 2 bin centers to construct edges.")
    edges = np.empty(c.size + 1, dtype=float)
    edges[1:-1] = 0.5 * (c[:-1] + c[1:])
    edges[0]    = c[0] - 0.5 * (c[1] - c[0])
    edges[-1]   = c[-1] + 0.5 * (c[-1] - c[-2])
    return edges


# -----------------------------
# Two-state projection (open / closed)
# -----------------------------

def compute_two_state_from_data(
    x_by_window,
    hb_by_window,
    bin_centers_sim,
    P_r,
    xmin_nm=None,
    xmax_nm=None,
    open_hb=OPEN_HB,
    closed_hb=CLOSED_HB,
    out_dat=OUT_OPEN_CLOSED_DAT,
):
    """
    Project WHAM P(r) onto two states defined by HB counts:

      open   : HB == open_hb
      closed : HB == closed_hb

    Uses existing WHAM P(r) (no WHAM recomputation).

    Distance restriction:
      - If xmin_nm / xmax_nm are provided, only bins with x_nm in [xmin_nm, xmax_nm]
        are considered, and frames whose distance falls outside this range are ignored.

    Returns a dict with:
      P_open, P_closed, dF_kBT, arrays N_tot, N_open, N_closed,
      P_open_given_r, P_closed_given_r, and mask 'use_bin' for used bins.
    """
    bin_centers_sim = np.asarray(bin_centers_sim)
    P_r             = np.asarray(P_r)

    # Bin edges and widths in simulation units
    edges_sim = construct_bin_edges_from_centers(bin_centers_sim)
    dx        = np.diff(edges_sim)  # length = nbins

    nbins = bin_centers_sim.size

    # Distance window in sim units (if provided)
    if xmin_nm is not None:
        xmin_sim = xmin_nm / SIM_TO_NM
    else:
        xmin_sim = edges_sim[0]

    if xmax_nm is not None:
        xmax_sim = xmax_nm / SIM_TO_NM
    else:
        xmax_sim = edges_sim[-1]

    # Mask of bins whose centers are within [xmin_sim, xmax_sim]
    use_bin = (bin_centers_sim >= xmin_sim) & (bin_centers_sim <= xmax_sim)
    if not np.any(use_bin):
        raise ValueError(
            f"No WHAM bins found in the requested range: "
            f"[{xmin_nm} nm, {xmax_nm} nm]."
        )

    N_tot    = np.zeros(nbins, dtype=np.int64)
    N_open   = np.zeros(nbins, dtype=np.int64)
    N_closed = np.zeros(nbins, dtype=np.int64)

    # Loop over windows and bin per-frame distances + HB counts
    for distances, hb_arr in zip(x_by_window, hb_by_window):
        if hb_arr is None:
            continue

        n = min(len(distances), len(hb_arr))
        if n == 0:
            continue

        d = np.asarray(distances[:n])
        h = np.asarray(hb_arr[:n])

        # Only keep frames whose distance is inside [xmin_sim, xmax_sim]
        mask_range = (d >= xmin_sim) & (d <= xmax_sim)
        if not np.any(mask_range):
            continue

        d = d[mask_range]
        h = h[mask_range]

        # Assign each frame to a bin index
        bin_indices = np.searchsorted(edges_sim, d, side="right") - 1
        mask_bins   = (bin_indices >= 0) & (bin_indices < nbins)
        if not np.any(mask_bins):
            continue

        bidx = bin_indices[mask_bins]
        h    = h[mask_bins]

        # Only collect into bins within use_bin
        valid_bins_mask = use_bin[bidx]
        if not np.any(valid_bins_mask):
            continue

        bidx = bidx[valid_bins_mask]
        h    = h[valid_bins_mask]

        # Total frames per bin
        N_tot += np.bincount(bidx, minlength=nbins)

        # Open frames per bin
        open_mask = (h == open_hb)
        if np.any(open_mask):
            N_open += np.bincount(bidx[open_mask], minlength=nbins)

        # Closed frames per bin
        closed_mask = (h == closed_hb)
        if np.any(closed_mask):
            N_closed += np.bincount(bidx[closed_mask], minlength=nbins)

    # Conditional probabilities P(open|r), P(closed|r)
    P_open_given_r   = np.zeros(nbins, dtype=float)
    P_closed_given_r = np.zeros(nbins, dtype=float)

    nonzero = (N_tot > 0) & use_bin
    P_open_given_r[nonzero]   = N_open[nonzero]   / N_tot[nonzero]
    P_closed_given_r[nonzero] = N_closed[nonzero] / N_tot[nonzero]

    # Use only bins in the selected range for integration
    P_r_eff = P_r.copy()
    P_r_eff[~use_bin] = 0.0

    # For integration, associate dx[j] = edges[j+1] - edges[j] with bin j
    dx_per_bin = dx  # length nbins, same as P_r and P_open_given_r

    P_open_tilde   = np.sum(P_r_eff * P_open_given_r   * dx_per_bin)
    P_closed_tilde = np.sum(P_r_eff * P_closed_given_r * dx_per_bin)
    P_sum          = P_open_tilde + P_closed_tilde

    if P_sum <= 0.0:
        raise RuntimeError(
            "No open or closed frames found in any bin within the chosen range."
        )

    P_open   = P_open_tilde   / P_sum
    P_closed = P_closed_tilde / P_sum

    dF_kBT = -np.log(P_closed / P_open)

    # Save per-bin profile (only bins in the selected range)
    if out_dat is not None:
        with open(out_dat, "w") as fh:
            fh.write("# x_sim  x_nm  P_r  P_open_given_r  P_closed_given_r  "
                     "N_tot  N_open  N_closed\n")
            for xs, Pr, po, pc, nt, no, nc, use in zip(
                bin_centers_sim,
                P_r,
                P_open_given_r,
                P_closed_given_r,
                N_tot,
                N_open,
                N_closed,
                use_bin,
            ):
                if not use:
                    continue
                fh.write(
                    f"{xs:.6f}  {xs*SIM_TO_NM:.6f}  {Pr:.6e}  "
                    f"{po:.6f}  {pc:.6f}  {nt:d}  {no:d}  {nc:d}\n"
                )
        print(f"[two-state] Wrote {out_dat}")

    print("\n[two-state] === Open/closed populations "
          f"(HB = {open_hb} vs {closed_hb}) ===")
    print(f"  P_open   = {P_open:.4f}")
    print(f"  P_closed = {P_closed:.4f}")
    print(f"  ΔF = F_closed - F_open = {dF_kBT:.3f} k_B T")
    print(f"      ≈ {dF_kBT * RT_KJMOL:.3f} kJ/mol  (T = {TEMPERATURE_K} K)")
    print("  (Populations are normalized within the {open, closed} subspace.)")
    print("=========================================\n")

    return {
        "P_open": P_open,
        "P_closed": P_closed,
        "dF_kBT": dF_kBT,
        "N_tot": N_tot,
        "N_open": N_open,
        "N_closed": N_closed,
        "P_open_given_r": P_open_given_r,
        "P_closed_given_r": P_closed_given_r,
        "use_bin": use_bin,
    }


# -----------------------------
# Plotting
# -----------------------------

def plot_two_state_profile(
    x_nm,
    F,
    P_open_given_r,
    P_closed_given_r,
    use_bin,
    out_png=OUT_OPEN_CLOSED_PNG,
):
    """
    Plot PMF (F) and P(open|r), P(closed|r) vs distance for bins in use_bin.
    """
    x_nm = np.asarray(x_nm)
    F    = np.asarray(F)

    mask = use_bin & np.isfinite(F)

    if not np.any(mask):
        print("[plot_two_state_profile] No valid bins to plot.")
        return

    x_nm_plot = x_nm[mask]
    F_plot    = F[mask]
    Po_plot   = P_open_given_r[mask]
    Pc_plot   = P_closed_given_r[mask]

    fig, ax1 = plt.subplots(figsize=(9, 5))

    color_pmf = "steelblue"
    ax1.plot(x_nm_plot, F_plot, "-", color=color_pmf, linewidth=2, label="PMF (WHAM)")
    ax1.set_xlabel("Anchor–anchor distance (nm)")
    ax1.set_ylabel("Free energy F (k_B T, from WHAM)", color=color_pmf)
    ax1.tick_params(axis="y", labelcolor=color_pmf)

    ax2 = ax1.twinx()
    color_open   = "tab:green"
    color_closed = "tab:red"

    ax2.plot(x_nm_plot, Po_plot, "--", color=color_open,  linewidth=2, label="P(open | r)")
    ax2.plot(x_nm_plot, Pc_plot, "-",  color=color_closed, linewidth=2, label="P(closed | r)")
    ax2.set_ylabel("State probability", color="black")
    ax2.tick_params(axis="y", labelcolor="black")
    ax2.set_ylim(-0.05, 1.05)

    l1, lb1 = ax1.get_legend_handles_labels()
    l2, lb2 = ax2.get_legend_handles_labels()
    ax1.legend(l1 + l2, lb1 + lb2, loc="upper right")

    plt.title("Two-state projection: PMF and P(open/closed | r)")
    fig.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()
    print(f"[plot] Saved {out_png}")


# -----------------------------
# Main
# -----------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compute two-state (open/closed) free energy "
                    "from existing WHAM profile + HB counts."
    )
    parser.add_argument(
        "--wham-file",
        required=True,
        help="Existing WHAM profile file (columns: x_sim x_nm F(kBT) P).",
    )
    parser.add_argument(
        "--xmin-nm",
        type=float,
        default=None,
        help="Minimum distance (nm) to include in analysis (default: WHAM min).",
    )
    parser.add_argument(
        "--xmax-nm",
        type=float,
        default=None,
        help="Maximum distance (nm) to include in analysis (default: WHAM max).",
    )
    args = parser.parse_args()

    print("=== 2state.py ===")
    print(f"CWD         : {os.getcwd()}")
    print(f"Temperature : {TEMPERATURE_K} K  (RT = {RT_KJMOL:.3f} kJ/mol)")
    print(f"WHAM file   : {args.wham_file}")
    print(f"xmin (nm)   : {args.xmin_nm}")
    print(f"xmax (nm)   : {args.xmax_nm}")

    umbrella_root = os.path.abspath(os.getcwd())

    # Load WHAM profile (no recomputation)
    x_sim_wham, x_nm_wham, F_wham, P_wham = load_wham_profile(args.wham_file)

    # Load KL indices for HB counting
    kl_indices = load_kl_indices(KL_CSV)

    # Find completed windows
    completed_windows = find_completed_windows(WINDOW_PATTERN, PROD_SUBDIR, TRAJ_FILENAME)
    if not completed_windows:
        print("ERROR: No completed windows found.")
        return

    print(f"\n[main] {len(completed_windows)} completed window(s).")

    x_by_window  = []
    hb_by_window = []

    # Collect distances and HB counts per window
    for w in completed_windows:
        prod_dir    = os.path.abspath(os.path.join(w, PROD_SUBDIR))
        traj_file   = os.path.join(prod_dir, TRAJ_FILENAME)
        forces_file = os.path.join(prod_dir, FORCES_FILENAME)

        print(f"\n[main] ── {w} ──")

        # Need anchors to compute distances
        r0, stiff, anchor_a, anchor_b = get_r0_stiff_anchors(forces_file)

        # Distances from trajectory
        distances = []
        for frame_i, positions in enumerate(
                parse_oxdna_traj_stream(traj_file, FRAME_STRIDE)):
            if frame_i >= MAX_FRAMES_PER_WINDOW:
                print(f"  Capped at {MAX_FRAMES_PER_WINDOW} frames.")
                break
            N = positions.shape[0]
            if anchor_a >= N or anchor_b >= N:
                raise IndexError(
                    f"Anchor index out of range (N={N}, "
                    f"anchor_a={anchor_a}, anchor_b={anchor_b})"
                )
            distances.append(float(np.linalg.norm(
                positions[anchor_a] - positions[anchor_b])))
            if (frame_i + 1) % 500 == 0:
                print(f"  {frame_i+1} distance frames processed...")

        distances = np.array(distances)
        if distances.size == 0:
            print(f"  WARNING: no distance frames parsed, skipping {w}.")
            continue

        x_by_window.append(distances)
        print(f"  Dist frames  : {len(distances)}")
        print(f"  Dist range   : {distances.min():.4f}–{distances.max():.4f} "
              f"(sim units)  mean={distances.mean():.4f} "
              f"({distances.mean()*SIM_TO_NM:.4f} nm)")

        # HB via OAT hb_list
        print(f"  [HB-OAT] Running hb_list for {w} ...")
        hb_arr = compute_hb_oat(prod_dir, INPUT_FILENAME, TRAJ_FILENAME,
                                kl_indices, ncpus=HB_NCPUS)
        hb_by_window.append(hb_arr)

    if not x_by_window:
        print("ERROR: No distance/HB data collected.")
        return

    # Make sure we write outputs to the umbrella root
    os.chdir(umbrella_root)

    # Two-state projection
    results = compute_two_state_from_data(
        x_by_window=x_by_window,
        hb_by_window=hb_by_window,
        bin_centers_sim=x_sim_wham,
        P_r=P_wham,
        xmin_nm=args.xmin_nm,
        xmax_nm=args.xmax_nm,
        open_hb=OPEN_HB,
        closed_hb=CLOSED_HB,
        out_dat=OUT_OPEN_CLOSED_DAT,
    )

    # Plot
    plot_two_state_profile(
        x_nm=x_nm_wham,
        F=F_wham,
        P_open_given_r=results["P_open_given_r"],
        P_closed_given_r=results["P_closed_given_r"],
        use_bin=results["use_bin"],
        out_png=OUT_OPEN_CLOSED_PNG,
    )

    print("\n[main] Done.")
    print(f"  {OUT_OPEN_CLOSED_DAT}  |  {OUT_OPEN_CLOSED_PNG}")


if __name__ == "__main__":
    main()
