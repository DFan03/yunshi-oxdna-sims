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
# Bend + twist from four segments
# ------------------------------------------------------------

def bend_twist_from_segments(LH1_coords, SH1_coords, SH2_coords, LH4_coords):
    """
    Compute bend and twist from segment coordinates.

    LH1_coords, SH1_coords, SH2_coords, LH4_coords:
        arrays of shape (Ni, 3) with Ni >= 2 for each segment.

    Bend:
      - Use LH1 and LH4 only (long helix).
      - Fit axis_LH1 and axis_LH4 using all specified residues (both strands).
      - bend_deg = angle between axis_LH1 and axis_LH4.

    Twist:
      - Treat each short helix as a dot (its centroid).
      - Treat long helix as a line: axis_L through C_L (using LH1+LH4).
      - Compute radial vectors from C_L to each dot, projected perpendicular to axis_L.
      - twist_deg = angle between those two radial directions (0–180°).
    """
    LH1_coords = np.asarray(LH1_coords)
    LH4_coords = np.asarray(LH4_coords)
    SH1_coords = np.asarray(SH1_coords)
    SH2_coords = np.asarray(SH2_coords)

    # ----- Bend: long helix only (LH1 vs LH4) -----
    axis_LH1 = fit_axis(LH1_coords)
    axis_LH4 = fit_axis(LH4_coords)
    bend_deg = angle_between(axis_LH1, axis_LH4)

    # ----- Long‑helix axis and centroid (LH1 + LH4 together) -----
    L_all_coords = np.vstack([LH1_coords, LH4_coords])
    axis_L = fit_axis(L_all_coords)
    axis_L /= np.linalg.norm(axis_L)
    C_L = L_all_coords.mean(axis=0)      # central point on the long helix

    # ----- Short‑helix centroids = the "dots" -----
    C1 = SH1_coords.mean(axis=0)
    C2 = SH2_coords.mean(axis=0)

    # ----- Radial vectors from long helix to each short helix -----
    R1 = C1 - C_L
    R2 = C2 - C_L

    # Remove components along the long‑helix axis → keep only perpendicular parts
    R1_perp = R1 - np.dot(R1, axis_L) * axis_L
    R2_perp = R2 - np.dot(R2, axis_L) * axis_L

    n1 = np.linalg.norm(R1_perp)
    n2 = np.linalg.norm(R2_perp)

    if n1 < 1e-8 or n2 < 1e-8:
        # One dot is almost on the axis; twist around the axis is undefined
        twist_deg = np.nan
    else:
        R1_perp /= n1
        R2_perp /= n2
        c = np.dot(R1_perp, R2_perp)
        c = np.clip(c, -1.0, 1.0)
        twist_deg = np.degrees(np.arccos(c))   # 0–180°

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
                    # Completed first frame
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
# Main (multimer-aware)
# ------------------------------------------------------------

def main():
    # Now: traj.dat junctions.txt k [exclude_fl]
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

    # Expect these four keys for ONE monomer definition
    required = ["LH1", "SH1", "SH2", "LH4"]
    for key in required:
        if key not in junctions_base:
            print(f"Error: junctions.txt must contain a line for '{key}'")
            sys.exit(1)

    # Determine total number of nucleotides and per-monomer length
    total_nt = get_total_nt(traj_file)
    if total_nt % k != 0:
        print(f"Warning: total_nt={total_nt} is not divisible by k={k}. "
              f"Using integer division for monomer length.")
    monomer_len = total_nt // k

    # Build segment indices for each monomer: shift by m * monomer_len
    # junctions_monomers[m][segment_name] = list of indices for monomer m
    junctions_monomers = []
    for m in range(k):
        offset = m * monomer_len
        jm = {}
        for name, base_indices in junctions_base.items():
            jm[name] = [i + offset for i in base_indices]
        junctions_monomers.append(jm)

    # Collect all nucleotide indices used across all monomers and segments
    all_indices = sorted({
        idx
        for jm in junctions_monomers
        for indices in jm.values()
        for idx in indices
    })

    # Open per-monomer CSVs
    monomer_files = []
    monomer_writers = []
    for m in range(k):
        fname = f"monomer_{m+1}.csv"
        f = open(fname, "w", newline="")
        w = csv.writer(f)
        w.writerow(["t", "twist_deg", "bend_deg"])
        monomer_files.append(f)
        monomer_writers.append(w)

    # Average CSV
    avg_file = open("average.csv", "w", newline="")
    avg_writer = csv.writer(avg_file)
    avg_writer.writerow(["t", "twist_avg_deg", "bend_avg_deg"])

    # Process trajectory
    for t, coords_all in parse_traj(traj_file, all_indices):
        coords_dict = {idx: coords for idx, coords in zip(all_indices, coords_all)}

        twists = []
        bends = []

        # Equivalent to running the original script once per monomer
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

            monomer_writers[m].writerow([t, twist, bend])
            twists.append(twist)
            bends.append(bend)

        # Average across monomers
        if exclude_fl and k > 2:
            use_twists = twists[1:-1]
            use_bends  = bends[1:-1]
        else:
            use_twists = twists
            use_bends  = bends

        if len(use_twists) > 0:
            twist_avg = float(np.nanmean(use_twists))
            bend_avg  = float(np.nanmean(use_bends))
            avg_writer.writerow([t, twist_avg, bend_avg])
        # If k<=2 and exclude_fl=True, no average is written for that timestep.

    # Close files
    for f in monomer_files:
        f.close()
    avg_file.close()

    print("Wrote monomer_*.csv and average.csv")

if __name__ == "__main__":
    main()
