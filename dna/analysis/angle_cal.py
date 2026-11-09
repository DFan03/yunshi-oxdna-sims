"""
oxDNA trajectory angle calculator
----------------------------------
Usage:
    python angle_pca.py trj.dat helices.csv

Output:
    angles.csv  —  columns: timestep, angle_deg

helices.csv format (strictly two lines):
    astart,aend
    bstart,bend
The order of start,end defines the intended direction of the helix vector:
the PCA axis is flipped if needed so that it points from the start-end
nucleotide toward the end nucleotide.

oxDNA trajectory frame format:
    t = <timestep>
    b = ...         (box, ignored)
    E = ...         (energy, ignored)
    <one line per nucleotide: rx ry rz ...>
"""

import sys
import csv
import numpy as np


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_helices(csv_path):
    """
    Return two (start, end, indices) tuples from a two-line start,end CSV.
    indices is a list of all nt indices in the range (inclusive),
    regardless of order, for coordinate extraction.
    start and end are kept separately for direction reference.
    """
    helices = []
    with open(csv_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(",")
            if len(parts) != 2:
                raise ValueError(
                    f"Each line in helices.csv must have exactly two values "
                    f"(start, end). Got: '{line}'"
                )
            start, end = int(parts[0]), int(parts[1])
            # Build index list in ascending order for array slicing;
            # direction is handled separately via start/end positions.
            lo, hi = min(start, end), max(start, end)
            indices = list(range(lo, hi + 1))
            helices.append((start, end, indices))
    if len(helices) != 2:
        raise ValueError(
            f"helices.csv must have exactly 2 non-empty lines, got {len(helices)}."
        )
    return helices[0], helices[1]


def pca_axis(coords, pos_start, pos_end):
    """
    Given an (N, 3) array of coordinates, return the first principal
    component (unit vector along the long axis of the helix), signed so
    that it points from pos_start toward pos_end.

    pos_start and pos_end are the 3D coordinates of the start and end
    nucleotides respectively (used only for sign disambiguation).
    """
    centred = coords - coords.mean(axis=0)
    _, _, Vt = np.linalg.svd(centred, full_matrices=False)
    axis = Vt[0]  # unit vector, sign arbitrary

    # Sign convention: axis should point in the same general direction
    # as the vector from start nucleotide to end nucleotide.
    reference = pos_end - pos_start
    if np.dot(axis, reference) < 0:
        axis = -axis
    return axis


def angle_between(v1, v2):
    """
    Return the angle in degrees between two vectors.
    Always returns a value in [0, 180].
    """
    cos_theta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    return np.degrees(np.arccos(cos_theta))


# ---------------------------------------------------------------------------
# oxDNA trajectory parser
# ---------------------------------------------------------------------------

def iter_frames(trj_path):
    """
    Generator that yields (timestep, positions) for each frame.

    Each frame starts with three header lines:
        t = <timestep>
        b = ...
        E = ...
    Followed by one line per nucleotide; first three values are rx, ry, rz.
    """
    with open(trj_path, "r") as f:
        lines = f.readlines()

    i = 0
    n_lines = len(lines)

    while i < n_lines:
        line = lines[i].strip()

        if not line.startswith("t ="):
            i += 1
            continue

        timestep = int(line.split("=")[1].strip())
        i += 1  # past 't = ...'
        i += 2  # skip 'b = ...' and 'E = ...'

        if i >= n_lines:
            break

        positions = []
        while i < n_lines:
            nt_line = lines[i].strip()
            if nt_line.startswith("t ="):
                break
            if nt_line:
                vals = nt_line.split()
                positions.append([float(vals[0]), float(vals[1]), float(vals[2])])
            i += 1

        if positions:
            yield timestep, np.array(positions)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) != 3:
        print("Usage: python angle_pca.py trj.dat helices.csv")
        sys.exit(1)

    trj_path = sys.argv[1]
    csv_path = sys.argv[2]

    h1_start, h1_end, helix1_idx = parse_helices(csv_path)[0]
    h2_start, h2_end, helix2_idx = parse_helices(csv_path)[1]

    print(f"Helix 1: nucleotides {h1_start} -> {h1_end} ({len(helix1_idx)} nts)")
    print(f"Helix 2: nucleotides {h2_start} -> {h2_end} ({len(helix2_idx)} nts)")

    results = []

    for frame_num, (timestep, positions) in enumerate(iter_frames(trj_path)):
        n_nt = len(positions)

        for idx in helix1_idx + helix2_idx:
            if idx >= n_nt:
                raise IndexError(
                    f"Frame {frame_num} (t={timestep}): index {idx} out of range "
                    f"(trajectory has {n_nt} nucleotides)."
                )

        coords1 = positions[helix1_idx]
        coords2 = positions[helix2_idx]

        # Direction reference: position of start nt and end nt in this frame
        vec1 = pca_axis(coords1, positions[h1_start], positions[h1_end])
        vec2 = pca_axis(coords2, positions[h2_start], positions[h2_end])

        theta = angle_between(vec1, vec2)
        results.append((timestep, theta))

    out_path = "angles.csv"
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["timestep", "angle_deg"])
        writer.writerows(results)

    print(f"Done. Wrote {len(results)} frames to {out_path}.")


if __name__ == "__main__":
    main()

