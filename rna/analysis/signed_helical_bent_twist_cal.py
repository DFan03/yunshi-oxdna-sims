#!/usr/bin/env python3
import sys
import numpy as np
import csv

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

    Steps:
      1) Subtract centroid.
      2) SVD of the centered coordinates.
      3) First row of Vt is principal direction.
      4) Normalize to length 1.
    """
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

def angle_between(u, v):
    """
    Angle in degrees between two vectors u, v (0–180°).
    """
    u = np.asarray(u)
    v = np.asarray(v)
    nu = np.linalg.norm(u)
    nv = np.linalg.norm(v)
    if nu < 1e-8 or nv < 1e-8:
        return np.nan
    u = u / nu
    v = v / nv
    c = np.dot(u, v)
    c = np.clip(c, -1.0, 1.0)
    return np.degrees(np.arccos(c))

# ------------------------------------------------------------
# Global bend reference normal
# ------------------------------------------------------------

def compute_global_ref_normal(LH1_coords, SH1_coords, SH2_coords):
    """
    Compute a global reference normal vector for the bend sign.

    We use the same construction as before (LH1 axis + SH1/SH2) but only ONCE,
    then keep this normal fixed for all frames.

    Returns:
        ref_normal (3,) unit vector
    """
    LH1_coords = np.asarray(LH1_coords)
    SH1_coords = np.asarray(SH1_coords)
    SH2_coords = np.asarray(SH2_coords)

    axis_LH1 = fit_axis(LH1_coords)

    C1 = SH1_coords.mean(axis=0)
    C2 = SH2_coords.mean(axis=0)
    diff = C2 - C1

    # Project diff perpendicular to LH1 axis
    diff_perp = diff - np.dot(diff, axis_LH1) * axis_LH1

    if np.linalg.norm(diff_perp) < 1e-8:
        # Degenerate: fall back to a fixed global axis (e.g., z-axis)
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

# ------------------------------------------------------------
# Bend + twist from four segments with GLOBAL bend sign
# ------------------------------------------------------------

def bend_twist_from_segments(LH1_coords, SH1_coords, SH2_coords, LH4_coords,
                             ref_normal_global):
    """
    Compute signed bend and twist from segment coordinates.

    Returns values in the range −180° to +180°.

    Twist sign convention:
        Positive = R1_perp → R2_perp is a right-hand rotation around axis_L.
        Sign determined by: sign(axis_L · (R1_perp × R2_perp))

    Bend magnitude:
        Angle between axis_LH1 and axis_LH4 (0–180°).

    Bend sign (GLOBAL):
        Uses a single reference normal vector 'ref_normal_global' that is
        fixed for the entire analysis. For each frame we compute:

            cross_axes = axis_LH1 × axis_LH4

        and set sign according to:

            if ref_normal_global · cross_axes < 0: bend = -bend

        This keeps the sign consistent across frames and avoids the
        per-frame sign determined by SH1/SH2.
    """
    LH1_coords = np.asarray(LH1_coords)
    LH4_coords = np.asarray(LH4_coords)
    SH1_coords = np.asarray(SH1_coords)
    SH2_coords = np.asarray(SH2_coords)

    # ----- Bend: long helix only (LH1 vs LH4) -----
    axis_LH1 = fit_axis(LH1_coords)
    axis_LH4 = fit_axis(LH4_coords)

    # Orient both axes so they roughly point from LH1 centroid toward LH4 centroid
    reference_dir = LH4_coords.mean(axis=0) - LH1_coords.mean(axis=0)
    if np.dot(axis_LH1, reference_dir) < 0:
        axis_LH1 = -axis_LH1
    if np.dot(axis_LH4, reference_dir) < 0:
        axis_LH4 = -axis_LH4

    # Unsigned bend angle (0–180)
    cos_bend = np.clip(np.dot(axis_LH1, axis_LH4), -1.0, 1.0)
    bend_deg = np.degrees(np.arccos(cos_bend))

    # Sign of bend: GLOBAL reference normal
    cross_axes = np.cross(axis_LH1, axis_LH4)
    if np.linalg.norm(cross_axes) > 1e-8:
        if np.dot(ref_normal_global, cross_axes) < 0:
            bend_deg = -bend_deg
    # If cross_axes is ~0 (perfectly straight), sign doesn't matter

    # ----- Long-helix axis and centroid (LH1 + LH4 together) -----
    L_all_coords = np.vstack([LH1_coords, LH4_coords])
    axis_L = fit_axis(L_all_coords)
    axis_L /= np.linalg.norm(axis_L)
    if np.dot(axis_L, axis_LH1) < 0:
        axis_L = -axis_L

    C_L = L_all_coords.mean(axis=0)

    # ----- Radial vectors from long helix to each short helix -----
    C1 = SH1_coords.mean(axis=0)
    C2 = SH2_coords.mean(axis=0)

    R1 = C1 - C_L
    R2 = C2 - C_L

    # Remove components along the long-helix axis → keep only perpendicular parts
    R1_perp = R1 - np.dot(R1, axis_L) * axis_L
    R2_perp = R2 - np.dot(R2, axis_L) * axis_L

    n1 = np.linalg.norm(R1_perp)
    n2 = np.linalg.norm(R2_perp)

    if n1 < 1e-8 or n2 < 1e-8:
        twist_deg = np.nan
    else:
        R1_perp /= n1
        R2_perp /= n2

        # Unsigned twist angle
        c = np.clip(np.dot(R1_perp, R2_perp), -1.0, 1.0)
        twist_deg = np.degrees(np.arccos(c))

        # Sign: positive = right-hand rotation of R1_perp → R2_perp around axis_L
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
    coords_dict = {}  # maps nt index -> np.array
    current_t = None

    with open(traj_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("t ="):
                # yield previous timestep if available
                if current_t is not None:
                    coords_list = [coords_dict[i] for i in coords_needed]
                    yield current_t, coords_list
                # new timestep
                current_t = float(line.split("=")[1].strip())
                coords_dict = {}
                nt_counter = 0
            elif line == "" or line.startswith("b") or line.startswith("E"):
                continue
            else:
                # parse coordinates
                parts = line.split()
                xyz = np.array([float(parts[0]), float(parts[1]), float(parts[2])])
                coords_dict[nt_counter] = xyz
                nt_counter += 1
        # yield last timestep
        if current_t is not None:
            coords_list = [coords_dict[i] for i in coords_needed]
            yield current_t, coords_list

# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} traj.dat junctions.txt")
        sys.exit(1)

    traj_file = sys.argv[1]
    junction_file = sys.argv[2]

    junctions = read_junctions(junction_file)

    # Expect these four keys
    required = ["LH1", "SH1", "SH2", "LH4"]
    for key in required:
        if key not in junctions:
            print(f"Error: junctions.txt must contain a line for '{key}'")
            sys.exit(1)

    # Collect all nucleotide indices used in any segment
    all_indices = sorted({i for indices in junctions.values() for i in indices})

    ref_normal_global = None  # will be set from the first frame

    # Output CSV
    with open("twist_bend.csv", "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["t", "twist_deg", "bend_deg"])

        for t, coords_all in parse_traj(traj_file, all_indices):
            # Map from global index -> coordinates for this frame
            coords_dict = {idx: coords for idx, coords in zip(all_indices, coords_all)}

            # Build per‑segment coordinate arrays in the specified index order
            LH1_coords = [coords_dict[i] for i in junctions["LH1"]]
            SH1_coords = [coords_dict[i] for i in junctions["SH1"]]
            SH2_coords = [coords_dict[i] for i in junctions["SH2"]]
            LH4_coords = [coords_dict[i] for i in junctions["LH4"]]

            # Initialize global bend reference normal from the FIRST frame
            if ref_normal_global is None:
                try:
                    ref_normal_global = compute_global_ref_normal(
                        LH1_coords, SH1_coords, SH2_coords
                    )
                except ValueError as e:
                    print(f"Error computing global bend reference normal: {e}")
                    sys.exit(1)

            try:
                twist, bend = bend_twist_from_segments(
                    LH1_coords, SH1_coords, SH2_coords, LH4_coords,
                    ref_normal_global
                )
            except ValueError:
                # Degenerate geometry for this frame
                twist, bend = np.nan, np.nan

            writer.writerow([t, twist, bend])

    print("Wrote twist_bend.csv")

if __name__ == "__main__":
    main()
