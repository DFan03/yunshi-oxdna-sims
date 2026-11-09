#!/usr/bin/env python3

import sys
import csv
import numpy as np

# ------------------------------------------------------------
# Read segment definitions
# ------------------------------------------------------------

def read_junctions(junction_file):
    """
    Read junctions.txt and return dictionary with nucleotide index lists (0-based).

    Expected format per line (no spaces), examples:

        LH1,1,10,101,110   # two ranges: 1–10 and 101–110
        SH1,20,30          # one range: 20–30
        SH2,50,60
        LH4,201,210,301,310

    Rules:
      - All indices are 1-based, inclusive.
      - If a line has 3 numbers (name,start,end)   → one range.
      - If a line has 5 numbers (name,s1,e1,s2,e2) → two ranges (concatenated).
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
                # one range
                start, end = nums
                if end < start:
                    raise ValueError(f"End < start in line: {line}")
                indices = list(range(start - 1, end))
            elif len(nums) == 4:
                # two ranges
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

# ------------------------------------------------------------
# Geometry helpers
# ------------------------------------------------------------

def fit_axis(points):
    """
    Given an array of shape (N, 3) with N >= 2, fit a best‑fit line and
    return a unit axis vector (directional vector) along that line.
    """
    P = np.asarray(points, dtype=float)
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

# ------------------------------------------------------------
# Local bend + twist from four segments (your original monomer geometry)
# ------------------------------------------------------------

def bend_twist_from_segments(LH1_coords, SH1_coords, SH2_coords, LH4_coords):
    """
    Compute local (per-monomer) twist and bend from the four segments:
    LH1, SH1, SH2, LH4.

    Returns:
        twist_deg, bend_deg (floats)
    """
    LH1_coords = np.asarray(LH1_coords, dtype=float)
    LH4_coords = np.asarray(LH4_coords, dtype=float)
    SH1_coords = np.asarray(SH1_coords, dtype=float)
    SH2_coords = np.asarray(SH2_coords, dtype=float)

    # ----- Bend: long helix only (LH1 vs LH4) -----
    axis_LH1 = fit_axis(LH1_coords)
    axis_LH4 = fit_axis(LH4_coords)

    # Orient axes consistently: both point from LH1 centroid toward LH4 centroid
    reference_dir = LH4_coords.mean(axis=0) - LH1_coords.mean(axis=0)
    if np.dot(axis_LH1, reference_dir) < 0:
        axis_LH1 = -axis_LH1
    if np.dot(axis_LH4, reference_dir) < 0:
        axis_LH4 = -axis_LH4

    # Unsigned bend, then sign it (local SH-based sign – we overwrite sign globally later)
    cos_bend = np.clip(np.dot(axis_LH1, axis_LH4), -1.0, 1.0)
    bend_deg = np.degrees(np.arccos(cos_bend))

    C1 = SH1_coords.mean(axis=0)
    C2 = SH2_coords.mean(axis=0)
    diff = C2 - C1
    diff_perp = diff - np.dot(diff, axis_LH1) * axis_LH1
    if np.linalg.norm(diff_perp) > 1e-8:
        ref_normal = np.cross(axis_LH1, diff_perp)
        ref_normal /= np.linalg.norm(ref_normal)
        cross_axes = np.cross(axis_LH1, axis_LH4)
        if np.dot(ref_normal, cross_axes) < 0:
            bend_deg = -bend_deg

    # ----- Long-helix axis and centroid (LH1 + LH4 together) -----
    L_all_coords = np.vstack([LH1_coords, LH4_coords])
    axis_L = fit_axis(L_all_coords)
    axis_L /= np.linalg.norm(axis_L)

    # Orient axis_L consistently with axis_LH1
    if np.dot(axis_L, axis_LH1) < 0:
        axis_L = -axis_L

    C_L = L_all_coords.mean(axis=0)

    # ----- Radial vectors from long helix to each short helix -----
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

# ------------------------------------------------------------
# Trajectory parser
# ------------------------------------------------------------

def parse_traj(traj_file, nt_indices):
    """
    Parse trajectory file and extract coordinates of the nucleotides in nt_indices.

    nt_indices: list of all nucleotide indices to extract (0-based).

    Yields:
        timestep (float),
        coords_list (list of np.array), in order of sorted(nt_indices).
    """
    coords_needed = sorted(nt_indices)
    coords_dict = {}
    current_t = None

    with open(traj_file) as f:
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

# ------------------------------------------------------------
# Helpers for multimer
# ------------------------------------------------------------

def get_total_nt(traj_file):
    """
    Count total number of nucleotides (coordinates) in one frame of traj_file.
    Assumes all frames have the same number of nucleotides.
    """
    with open(traj_file) as f:
        in_frame = False
        nt_counter = 0
        for line in f:
            line = line.strip()
            if line.startswith("t ="):
                if in_frame:
                    break
                else:
                    in_frame = True
                    nt_counter = 0
            elif not in_frame:
                continue
            elif line == "" or line.startswith("b") or line.startswith("E"):
                continue
            else:
                nt_counter += 1
        return nt_counter

def parse_bool(s):
    """
    Parse a string to bool. Accepts: true/false, 1/0, yes/no, y/n (case-insensitive)
    """
    s = s.lower()
    if s in ("1", "true", "yes", "y"):
        return True
    if s in ("0", "false", "no", "n"):
        return False
    raise ValueError(f"Cannot parse boolean value from '{s}'")

# ------------------------------------------------------------
# Backbone axis and CW/CCW sign around it
# ------------------------------------------------------------

def compute_monomer_centroids(coords_dict, junctions_monomers):
    """
    For one frame, compute centroid of LH1+LH4 for each monomer.

    Returns:
        C_L_array: ndarray shape (k, 3)
    """
    C_L_list = []
    for jm in junctions_monomers:
        LH1_coords = [coords_dict[i] for i in jm["LH1"]]
        LH4_coords = [coords_dict[i] for i in jm["LH4"]]
        L_all_coords = np.vstack([LH1_coords, LH4_coords])
        C_L = np.mean(L_all_coords, axis=0)
        C_L_list.append(C_L)
    return np.asarray(C_L_list, dtype=float)

def compute_axis_and_basis_from_first_last(C_L_array):
    """
    Define a global axis from first to last monomer and a local (e1, e2, axis_Z) basis.

    - axis_Z points from monomer 1 to monomer k.
    - e1 is a radial direction from axis to monomer 1.
    - e2 completes a right-handed basis.
    """
    C = np.asarray(C_L_array, dtype=float)
    k = C.shape[0]
    if k < 2:
        raise ValueError("Need at least 2 monomers to define axis")

    axis_Z = C[-1] - C[0]
    norm_Z = np.linalg.norm(axis_Z)
    if norm_Z < 1e-8:
        axis_Z = fit_axis(C)
    else:
        axis_Z /= norm_Z

    C_mean = C.mean(axis=0)

    R0 = C[0] - C_mean
    R0_perp = R0 - np.dot(R0, axis_Z) * axis_Z
    nR0 = np.linalg.norm(R0_perp)
    if nR0 < 1e-8 and k > 2:
        R1 = C[1] - C_mean
        R1_perp = R1 - np.dot(R1, axis_Z) * axis_Z
        nR1 = np.linalg.norm(R1_perp)
        if nR1 > 1e-8:
            R0_perp = R1_perp
            nR0 = nR1
    if nR0 < 1e-8:
        tmp = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(tmp, axis_Z)) > 0.9:
            tmp = np.array([0.0, 1.0, 0.0])
        R0_perp = tmp - np.dot(tmp, axis_Z) * axis_Z
        nR0 = np.linalg.norm(R0_perp)
        if nR0 < 1e-8:
            raise ValueError("Cannot define a radial direction for backbone basis")

    e1 = R0_perp / nR0
    e2 = np.cross(axis_Z, e1)
    e2 /= np.linalg.norm(e2)

    return axis_Z, C_mean, e1, e2

def compute_backbone_angles(C_L_array, axis_Z, C_mean, e1, e2):
    """
    Compute angular position theta_m (in radians) of each monomer around axis_Z,
    in the (e1, e2) plane.

    Returns:
        theta: array of shape (k,) with angles in (-pi, pi].
    """
    C = np.asarray(C_L_array, dtype=float)
    k = C.shape[0]
    theta = np.zeros(k)
    for m in range(k):
        R = C[m] - C_mean
        R_perp = R - np.dot(R, axis_Z) * axis_Z
        x = np.dot(R_perp, e1)
        y = np.dot(R_perp, e2)
        theta[m] = np.arctan2(y, x)
    return theta

def compute_cw_signs_from_theta(theta):
    """
    Given angular positions theta[m] (radians) along the chain, determine a CW/CCW
    sign for each monomer, with REVERSED convention:

      CW (dtheta < 0)  --> sign = -1
      CCW (dtheta > 0) --> sign = +1

    Implementation:

      - For each pair m -> m+1, compute dtheta (wrapped to (-pi, pi]).
      - For monomer indices:
          sign[0]      = sign_from_delta(dtheta[0])
          sign[k-1]    = sign_from_delta(dtheta[k-2])
          sign[m]      = sign_from_delta( (dtheta[m-1] + dtheta[m]) / 2 )  for 1<=m<=k-2

    Returns:
        sign: array of length k with values in {-1, 0, +1}.
    """
    theta = np.asarray(theta, dtype=float)
    k = theta.shape[0]
    sign = np.zeros(k)
    if k < 2:
        return sign

    dtheta = np.zeros(k - 1)
    for m in range(k - 1):
        d = theta[m + 1] - theta[m]
        while d > np.pi:
            d -= 2.0 * np.pi
        while d <= -np.pi:
            d += 2.0 * np.pi
        dtheta[m] = d

    def sign_from_delta(d):
        eps = 1e-3  # ~0.06 degrees
        if d < -eps:
            # CW step: dtheta < 0 --> NOW sign = -1
            return -1.0
        elif d > eps:
            # CCW step: dtheta > 0 --> NOW sign = +1
            return 1.0
        else:
            return 0.0

    sign[0] = sign_from_delta(dtheta[0])

    for m in range(1, k - 1):
        avg = 0.5 * (dtheta[m - 1] + dtheta[m])
        sign[m] = sign_from_delta(avg)

    sign[k - 1] = sign_from_delta(dtheta[k - 2])

    return sign

# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():
    # Usage: traj.dat junctions.txt k [exclude_fl]
    if len(sys.argv) < 4 or len(sys.argv) > 5:
        print(f"Usage: {sys.argv[0]} traj.dat junctions.txt k [exclude_fl]")
        print("  traj.dat       : trajectory file")
        print("  junctions.txt  : junction indices for a single monomer (1-based)")
        print("  k              : number of monomers (integer >= 1)")
        print("  exclude_fl     : optional, if 'true' or '1', exclude first and last monomer from average")
        sys.exit(1)

    traj_file = sys.argv[1]
    junction_file = sys.argv[2]

    try:
        k = int(sys.argv[3])
    except ValueError:
        print("Error: k must be an integer")
        sys.exit(1)

    if k < 1:
        print("Error: k must be >= 1")
        sys.exit(1)

    exclude_fl = False
    if len(sys.argv) == 5:
        try:
            exclude_fl = parse_bool(sys.argv[4])
        except ValueError as e:
            print(f"Error parsing exclude_fl: {e}")
            sys.exit(1)

    junctions_base = read_junctions(junction_file)

    required = ["LH1", "SH1", "SH2", "LH4"]
    for key in required:
        if key not in junctions_base:
            print(f"Error: junctions.txt must contain a line for '{key}'")
            sys.exit(1)

    total_nt = get_total_nt(traj_file)
    if total_nt % k != 0:
        print(f"Warning: total_nt={total_nt} is not divisible by k={k}. "
              f"Using integer division for monomer length.")
    monomer_len = total_nt // k

    junctions_monomers = []
    for m in range(k):
        offset = m * monomer_len
        jm = {}
        for name, base_indices in junctions_base.items():
            jm[name] = [i + offset for i in base_indices]
        junctions_monomers.append(jm)

    all_indices = sorted({
        idx
        for jm in junctions_monomers
        for indices in jm.values()
        for idx in indices
    })

    monomer_files = []
    monomer_writers = []
    for m in range(k):
        fname = f"monomer_{m+1}.csv"
        f = open(fname, "w", newline="")
        w = csv.writer(f)
        w.writerow(["t", "twist_deg", "bend_deg"])
        monomer_files.append(f)
        monomer_writers.append(w)

    avg_file = open("average.csv", "w", newline="")
    avg_writer = csv.writer(avg_file)
    avg_writer.writerow(["t", "twist_avg_deg", "bend_avg_deg"])

    # Process trajectory
    for t, coords_all in parse_traj(traj_file, all_indices):
        coords_dict = {idx: coords for idx, coords in zip(all_indices, coords_all)}

        # 1) Local twist/bend
        twist_local = [np.nan] * k
        bend_local = [np.nan] * k

        for m, jm in enumerate(junctions_monomers):
            LH1_coords = [coords_dict[i] for i in jm["LH1"]]
            SH1_coords = [coords_dict[i] for i in jm["SH1"]]
            SH2_coords = [coords_dict[i] for i in jm["SH2"]]
            LH4_coords = [coords_dict[i] for i in jm["LH4"]]

            try:
                twist, bend = bend_twist_from_segments(
                    LH1_coords, SH1_coords, SH2_coords, LH4_coords
                )
            except ValueError:
                twist, bend = np.nan, np.nan

            twist_local[m] = twist
            bend_local[m] = bend

        # 2) Global axis and CW/CCW sign
        C_L_array = compute_monomer_centroids(coords_dict, junctions_monomers)
        try:
            axis_Z, C_mean, e1, e2 = compute_axis_and_basis_from_first_last(C_L_array)
        except ValueError:
            twist_final = twist_local
            bend_final = [0.0 if np.isfinite(b) else np.nan for b in bend_local]
        else:
            theta = compute_backbone_angles(C_L_array, axis_Z, C_mean, e1, e2)
            sign_cw = compute_cw_signs_from_theta(theta)

            twist_final = twist_local
            bend_final = [np.nan] * k
            for m in range(k):
                s = sign_cw[m]
                if np.isfinite(s) and np.isfinite(bend_local[m]):
                    bend_final[m] = s * abs(bend_local[m])
                else:
                    bend_final[m] = np.nan

        # 3) Write per-monomer data
        for m in range(k):
            monomer_writers[m].writerow([t, twist_final[m], bend_final[m]])

        # 4) Averages
        if exclude_fl and k > 2:
            use_twist = twist_final[1:-1]
            use_bend = bend_final[1:-1]
        else:
            use_twist = twist_final
            use_bend = bend_final

        if len(use_twist) > 0:
            twist_avg = float(np.nanmean(use_twist))
            bend_avg = float(np.nanmean(use_bend))
            avg_writer.writerow([t, twist_avg, bend_avg])

    for f in monomer_files:
        f.close()
    avg_file.close()

    print("Wrote monomer_*.csv and average.csv with REVERSED CW/CCW-signed bend and original twist")

if __name__ == "__main__":
    main()
