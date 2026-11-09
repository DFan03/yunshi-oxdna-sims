#!/usr/bin/env python3
import sys
import numpy as np
import csv

def read_junctions(junction_file):
    """
    Read junctions.txt and return dictionary with nucleotide indices (0-based)
    """
    junctions = {}
    with open(junction_file) as f:
        for line in f:
            if line.strip() == "":
                continue
            parts = line.strip().split(",")
            name = parts[0].strip()
            indices = [int(x)-1 for x in parts[1:]]  # convert to 0-based
            junctions[name] = indices
    return junctions

def bend_twist_from_coords(coords):
    """
    coords: 8 points (np.array of shape (3,))
    order: a1,a2,b1,b2,c1,c2,d1,d2
    returns: bend, twist in degrees
    """
    a1,a2,b1,b2,c1,c2,d1,d2 = coords

    # Vectors along lines
    vA = a2 - a1
    vB = b2 - b1
    vC = c2 - c1
    vD = d2 - d1

    # Junction normals
    n1 = np.cross(vA, vB)
    n2 = np.cross(vC, vD)
    n1 /= np.linalg.norm(n1)
    n2 /= np.linalg.norm(n2)

    # Junction centers
    O1 = (a1 + a2 + b1 + b2)/4
    O2 = (c1 + c2 + d1 + d2)/4
    axis = O2 - O1
    axis /= np.linalg.norm(axis)

    # Bend angle
    bend_rad = np.arccos(np.clip(np.dot(n1,n2), -1.0, 1.0))
    bend_deg = np.degrees(bend_rad)

    # Project normals onto plane perpendicular to axis
    n1_proj = n1 - np.dot(n1, axis)*axis
    n2_proj = n2 - np.dot(n2, axis)*axis
    n1_proj /= np.linalg.norm(n1_proj)
    n2_proj /= np.linalg.norm(n2_proj)

    # Twist angle
    cos_phi = np.dot(n1_proj, n2_proj)
    sin_phi = np.dot(axis, np.cross(n1_proj, n2_proj))
    twist_rad = np.arctan2(sin_phi, cos_phi)
    twist_deg = np.degrees(twist_rad)

    return twist_deg, bend_deg

def parse_traj(traj_file, nt_indices):
    """
    Parse trajectory file and extract coordinates of the nucleotides in nt_indices
    nt_indices: list of all nucleotide indices to extract (0-based)
    yields: timestep, coords_list
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
                    # assemble coords in requested order
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

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} traj.dat junctions.txt")
        sys.exit(1)

    traj_file = sys.argv[1]
    junction_file = sys.argv[2]

    junctions = read_junctions(junction_file)

    # Order of nucleotides for bend_twist: a1,a2,b1,b2,c1,c2,d1,d2
    # So extract indices in that order
    nt_order = (
        junctions["a1"] +
        junctions["a2"] +
        junctions["b1"] +
        junctions["b2"]
    )

    # flatten to a list of unique nt indices
    nt_indices = sorted(set(nt_order))

    # output CSV
    with open("twist_bend.csv", "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["t", "twist_deg", "bend_deg"])

        for t, coords_all in parse_traj(traj_file, nt_indices):
            # map coords_all to the order expected by bend_twist_from_coords
            coords_dict = {idx: coords for idx, coords in zip(nt_indices, coords_all)}
            coords_ordered = [coords_dict[i] for i in nt_order]
            twist, bend = bend_twist_from_coords(coords_ordered)
            writer.writerow([t, twist, bend])

    print("Wrote twist_bend.csv")

if __name__ == "__main__":
    main()
