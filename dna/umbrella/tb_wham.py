#!/usr/bin/env python3
import sys
import os
import glob
import csv
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# User parameters
# ============================================================

T_K = 300.0                  # Temperature in K
kB_kcal = 0.0019872041       # kcal/mol/K
beta = 1.0 / (kB_kcal * T_K)

# Trajectories are in ./window_xx/prod/trajectory.dat
TRAJ_PATTERN = "window_*/prod/trajectory.dat"

# Binning for the apparent PMFs
NBINS_TWIST = 80
NBINS_BEND  = 80
MIN_TWIST = None   # if None → inferred from data
MAX_TWIST = None
MIN_BEND  = 0.0    # bend magnitude, so min=0 by default
MAX_BEND  = None

# ============================================================
# Read helix / segment definitions (helices.csv)
# ============================================================

def read_junctions(junction_file):
    """
    Read helices.csv and return dict name -> list of 0-based nucleotide indices.

    Expected format per line (CSV, no spaces), examples:

        LH1,1,10,101,110
        SH1,20,30
        SH2,50,60
        LH4,201,210,301,310
    """
    junctions = {}
    with open(junction_file) as f:
        for line in f:
            line = line.strip()
            if line == "":
                continue
            parts = [p.strip() for p in line.split(",")]
            name = parts[0]
            nums = [int(x) for x in parts[1:]]
            if len(nums) == 2:
                start, end = nums
                if end < start:
                    raise ValueError(f"End < start in line: {line}")
                indices = list(range(start - 1, end))
            elif len(nums) == 4:
                s1, e1, s2, e2 = nums
                if e1 < s1 or e2 < s2:
                    raise ValueError(f"End < start in line: {line}")
                indices1 = list(range(s1 - 1, e1))
                indices2 = list(range(s2 - 1, e2))
                indices = indices1 + indices2
            else:
                raise ValueError(
                    f"Line '{line}' must have 3 or 5 fields (name,start,end[,start2,end2])"
                )
            junctions[name] = indices
    return junctions

# ============================================================
# Geometry helpers (from your updated twist/bend implementation)
# ============================================================

def fit_axis(points):
    P = np.asarray(points)
    if P.shape[0] < 2:
        raise ValueError("Need at least 2 points to fit an axis")
    center = P.mean(axis=0)
    X = P - center
    _, _, Vt = np.linalg.svd(X)
    axis = Vt[0]
    norm = np.linalg.norm(axis)
    if norm < 1e-8:
        raise ValueError("Axis ill-defined (points nearly identical)")
    axis /= norm
    return axis

def compute_global_ref_normal(LH1_coords, SH1_coords, SH2_coords):
    LH1_coords = np.asarray(LH1_coords)
    SH1_coords = np.asarray(SH1_coords)
    SH2_coords = np.asarray(SH2_coords)

    axis_LH1 = fit_axis(LH1_coords)

    C1 = SH1_coords.mean(axis=0)
    C2 = SH2_coords.mean(axis=0)
    diff = C2 - C1

    diff_perp = diff - np.dot(diff, axis_LH1) * axis_LH1
    if np.linalg.norm(diff_perp) < 1e-8:
        ref_normal = np.array([0.0, 0.0, 1.0])
        ref_normal /= np.linalg.norm(ref_normal)
        print(
            "Warning: diff_perp nearly zero when computing global bend reference; "
            "falling back to global z-axis."
        )
        return ref_normal

    ref_normal = np.cross(axis_LH1, diff_perp)
    nrm = np.linalg.norm(ref_normal)
    if nrm < 1e-8:
        ref_normal = np.array([0.0, 0.0, 1.0])
        ref_normal /= np.linalg.norm(ref_normal)
        print(
            "Warning: reference normal ill-defined when computing global bend reference; "
            "falling back to global z-axis."
        )
        return ref_normal

    ref_normal /= nrm
    return ref_normal

def bend_twist_from_segments(LH1_coords, SH1_coords, SH2_coords, LH4_coords,
                             ref_normal_global):
    """
    Compute signed twist and signed bend from four segments.

    Returns:
        twist_deg ∈ (-180, 180]
        bend_deg  ∈ (-180, 180]  (we later use |bend_deg| for PMF)
    """
    LH1_coords = np.asarray(LH1_coords)
    LH4_coords = np.asarray(LH4_coords)
    SH1_coords = np.asarray(SH1_coords)
    SH2_coords = np.asarray(SH2_coords)

    # ----- Bend: long helix only (LH1 vs LH4) -----
    axis_LH1 = fit_axis(LH1_coords)
    axis_LH4 = fit_axis(LH4_coords)

    reference_dir = LH4_coords.mean(axis=0) - LH1_coords.mean(axis=0)
    if np.dot(axis_LH1, reference_dir) < 0:
        axis_LH1 = -axis_LH1
    if np.dot(axis_LH4, reference_dir) < 0:
        axis_LH4 = -axis_LH4

    cos_bend = np.clip(np.dot(axis_LH1, axis_LH4), -1.0, 1.0)
    bend_deg = np.degrees(np.arccos(cos_bend))

    cross_axes = np.cross(axis_LH1, axis_LH4)
    if np.linalg.norm(cross_axes) > 1e-8:
        if np.dot(ref_normal_global, cross_axes) < 0:
            bend_deg = -bend_deg

    # ----- Long-helix axis & centroid (LH1 + LH4) -----
    L_all_coords = np.vstack([LH1_coords, LH4_coords])
    axis_L = fit_axis(L_all_coords)
    axis_L /= np.linalg.norm(axis_L)
    if np.dot(axis_L, axis_LH1) < 0:
        axis_L = -axis_L
    C_L = L_all_coords.mean(axis=0)

    # ----- Short-helix centers -----
    C1 = SH1_coords.mean(axis=0)
    C2 = SH2_coords.mean(axis=0)

    R1 = C1 - C_L
    R2 = C2 - C_L

    R1_perp = R1 - np.dot(R1, axis_L) * axis_L
    R2_perp = R2 - np.dot(R2, axis_L) * axis_L

    n1 = np.linalg.norm(R1_perp)
    n2 = np.linalg.norm(R2_perp)

    if n1 < 1e-8 or n2 < 1e-8:
        twist_deg = np.nan
    else:
        R1_perp /= n1
        R2_perp /= n2

        c = np.clip(np.dot(R1_perp, R2_perp), -1.0, 1.0)
        twist_deg = np.degrees(np.arccos(c))

        cross = np.cross(R1_perp, R2_perp)
        if np.dot(axis_L, cross) < 0:
            twist_deg = -twist_deg

    return twist_deg, bend_deg

# ============================================================
# Trajectory parser
# ============================================================

def parse_traj(traj_file, nt_indices):
    """
    Parse trajectory file and extract coordinates of nt_indices (0-based).

    Expects:
      - lines starting with 't =' introduce new frame
      - coordinate lines: "x y z"
      - lines starting with 'b' or 'E' or empty are skipped.
    """
    coords_needed = sorted(nt_indices)
    coords_dict = {}
    current_t = None

    with open(traj_file) as f:
        nt_counter = 0
        for line in f:
            line = line.strip()
            if line.startswith("t ="):
                if current_t is not None:
                    coords_list = [coords_dict[i] for i in coords_needed]
                    yield current_t, coords_list
                current_t = float(line.split("=")[1].strip())
                coords_dict = {}
                nt_counter = 0
            elif line == "" or line.startswith("b") or line.startswith("E"):
                continue
            else:
                parts = line.split()
                xyz = np.array([float(parts[0]), float(parts[1]), float(parts[2])])
                coords_dict[nt_counter] = xyz
                nt_counter += 1

        if current_t is not None:
            coords_list = [coords_dict[i] for i in coords_needed]
            yield current_t, coords_list

# ============================================================
# Compute twist & bend series for one trajectory
# ============================================================

def compute_twist_bend_series(traj_file, junctions, ref_normal_global_container):
    all_indices = sorted({i for indices in junctions.values() for i in indices})
    twists = []
    bends  = []

    for t, coords_all in parse_traj(traj_file, all_indices):
        coords_dict = {idx: coords for idx, coords in zip(all_indices, coords_all)}

        LH1_coords = [coords_dict[i] for i in junctions["LH1"]]
        SH1_coords = [coords_dict[i] for i in junctions["SH1"]]
        SH2_coords = [coords_dict[i] for i in junctions["SH2"]]
        LH4_coords = [coords_dict[i] for i in junctions["LH4"]]

        if ref_normal_global_container["ref"] is None:
            try:
                ref = compute_global_ref_normal(LH1_coords, SH1_coords, SH2_coords)
            except ValueError as e:
                print(f"Error computing global bend reference normal: {e}")
                sys.exit(1)
            ref_normal_global_container["ref"] = ref

        try:
            twist, bend = bend_twist_from_segments(
                LH1_coords, SH1_coords, SH2_coords, LH4_coords,
                ref_normal_global_container["ref"]
            )
        except ValueError:
            twist, bend = np.nan, np.nan

        if not np.isnan(twist) and not np.isnan(bend):
            twists.append(twist)
            bends.append(bend)

    return np.array(twists), np.array(bends)

# ============================================================
# Plotting helpers
# ============================================================

def apparent_pmf_1d(values, nbins, vmin=None, vmax=None, label="CV",
                    out_dat="pmf.dat", out_png="pmf.png"):
    """
    F(x) = -k_B T ln P(x) (apparent, no reweighting by distance bias).
    """
    vals = np.asarray(values)
    if vmin is None:
        vmin = np.min(vals)
    if vmax is None:
        vmax = np.max(vals)

    hist, edges = np.histogram(vals, bins=nbins, range=(vmin, vmax))
    centers = 0.5 * (edges[:-1] + edges[1:])

    mask = hist > 0
    P = np.zeros_like(hist, dtype=float)
    P[mask] = hist[mask] / np.sum(hist[mask])

    F = np.full_like(P, np.nan, dtype=float)
    F[mask] = -kB_kcal * T_K * np.log(P[mask])
    if np.any(mask):
        F[mask] -= np.nanmin(F[mask])

    with open(out_dat, "w") as f:
        f.write(f"# {label}   F_kcal_per_mol   counts\n")
        for x, Fx, cnt in zip(centers, F, hist):
            Fx_out = "nan" if np.isnan(Fx) else f"{Fx:12.5f}"
            f.write(f"{x:12.5f}  {Fx_out}  {int(cnt)}\n")
    print(f"Wrote {out_dat}")

    plt.figure()
    mask_plot = ~np.isnan(F)
    plt.plot(centers[mask_plot], F[mask_plot], "-o", ms=3)
    plt.xlabel(label)
    plt.ylabel("Apparent free energy (kcal/mol)")
    plt.title(f"Apparent PMF along {label}")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()
    print(f"Wrote {out_png}")

def plot_violin_per_window(data_list, window_ids, centers_hint, ylabel, out_png, title):
    positions = np.arange(1, len(data_list) + 1)

    plt.figure(figsize=(8, 4))
    parts = plt.violinplot(
        data_list,
        positions=positions,
        showmeans=False,
        showmedians=True,
        showextrema=True,
    )

    for pc in parts["bodies"]:
        pc.set_facecolor("#1f77b4")
        pc.set_edgecolor("black")
        pc.set_alpha(0.6)
    for partname in ("cbars", "cmins", "cmaxes", "cmedians"):
        vp = parts[partname]
        vp.set_edgecolor("black")
        vp.set_linewidth(1)

    labels = [f"{wid} ({c:.1f})" for wid, c in zip(window_ids, centers_hint)]
    plt.xticks(positions, labels, rotation=45, ha="right")

    plt.xlabel("Window (id, mean value)")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()
    print(f"Wrote {out_png}")

# ============================================================
# Main
# ============================================================

def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} helices.csv")
        sys.exit(1)

    helices_file = sys.argv[1]

    junctions = read_junctions(helices_file)
    for key in ["LH1", "SH1", "SH2", "LH4"]:
        if key not in junctions:
            print(f"Error: {helices_file} must contain a line for '{key}'")
            sys.exit(1)

    # ./window_xx/prod/trajectory.dat
    traj_files = sorted(glob.glob(TRAJ_PATTERN))
    if len(traj_files) == 0:
        print(f"No trajectory files matching pattern '{TRAJ_PATTERN}'")
        sys.exit(1)

    # map window id -> trajectory
    traj_map = {}
    for fn in traj_files:
        # fn example: "window_00/prod/trajectory.dat"
        prod_dir   = os.path.dirname(fn)            # "window_00/prod"
        win_parent = os.path.dirname(prod_dir)      # "window_00"
        win_dir    = os.path.basename(win_parent)   # "window_00"

        if not win_dir.startswith("window_"):
            raise RuntimeError(
                f"Cannot parse window id from directory '{win_dir}' for file '{fn}'"
            )
        try:
            w = int(win_dir.split("window_")[1])
        except Exception:
            raise RuntimeError(
                f"Cannot parse integer window id from '{win_dir}' for file '{fn}'"
            )

        if w in traj_map:
            print(f"Warning: multiple trajectory.dat for window {w}; using the first one only.")
            continue
        traj_map[w] = fn

    window_ids_sorted = sorted(traj_map.keys())
    print(f"Found windows: {window_ids_sorted}")

    ref_normal_global_container = {"ref": None}
    twist_per_window = []
    bendmag_per_window = []

    for w in window_ids_sorted:
        fn = traj_map[w]
        print(f"Processing window {w}: {fn}")
        twist_arr, bend_arr = compute_twist_bend_series(
            fn, junctions, ref_normal_global_container
        )
        bend_mag = np.abs(bend_arr)

        print(f"  {len(twist_arr)} valid frames")

        twist_per_window.append(twist_arr)
        bendmag_per_window.append(bend_mag)

    all_twist = np.concatenate(twist_per_window)
    all_bend_mag = np.concatenate(bendmag_per_window)

    twist_min = MIN_TWIST if MIN_TWIST is not None else np.min(all_twist)
    twist_max = MAX_TWIST if MAX_TWIST is not None else np.max(all_twist)

    bend_min = MIN_BEND
    bend_max = MAX_BEND if MAX_BEND is not None else np.max(all_bend_mag)

    apparent_pmf_1d(
        all_twist,
        nbins=NBINS_TWIST,
        vmin=twist_min,
        vmax=twist_max,
        label="Twist (deg)",
        out_dat="pmf_twist.dat",
        out_png="pmf_twist.png",
    )

    apparent_pmf_1d(
        all_bend_mag,
        nbins=NBINS_BEND,
        vmin=bend_min,
        vmax=bend_max,
        label="|Bend| (deg)",
        out_dat="pmf_bend.dat",
        out_png="pmf_bend.png",
    )

    twist_means = [np.mean(tw) for tw in twist_per_window]
    bend_means  = [np.mean(bm) for bm in bendmag_per_window]

    plot_violin_per_window(
        twist_per_window,
        window_ids_sorted,
        twist_means,
        ylabel="Twist (deg)",
        out_png="twist_windows_violin.png",
        title="Twist distributions per window",
    )

    plot_violin_per_window(
        bendmag_per_window,
        window_ids_sorted,
        bend_means,
        ylabel="|Bend| (deg)",
        out_png="bend_windows_violin.png",
        title="Bend-magnitude distributions per window",
    )

if __name__ == "__main__":
    main()
