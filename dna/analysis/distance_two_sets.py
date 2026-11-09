import numpy as np
import csv
import sys

def parse_dat(filepath):
    """Parse .dat config file, extract XYZ coordinates from nucleotide lines."""
    coords = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip header lines (t, b, E)
            if not line or line.startswith('t') or line.startswith('b') or line.startswith('E'):
                continue
            values = line.split()
            if len(values) >= 3:
                try:
                    x, y, z = float(values[0]), float(values[1]), float(values[2])
                    coords.append([x, y, z])
                except ValueError:
                    continue
    return np.array(coords)

def parse_csv(filepath):
    """Parse CSV file with nucleotide indices, return as list of ints."""
    indices = []
    with open(filepath, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            for val in row:
                val = val.strip()
                if val:
                    indices.append(int(val))
    return indices

def compute_com(coords):
    """Compute center of mass (equal weights assumed)."""
    return np.mean(coords, axis=0)

def main():
    if len(sys.argv) != 4:
        print("Usage: python script.py a.csv b.csv x.dat")
        sys.exit(1)

    a_csv, b_csv, dat_file = sys.argv[1], sys.argv[2], sys.argv[3]

    # Load all nucleotide coordinates from .dat
    all_coords = parse_dat(dat_file)
    print(f"Total nucleotides in .dat: {len(all_coords)}")

    # Load indices from CSVs (0-based assumed)
    indices_a = parse_csv(a_csv)
    indices_b = parse_csv(b_csv)
    print(f"Nucleotides in group A: {len(indices_a)}")
    print(f"Nucleotides in group B: {len(indices_b)}")

    # Extract coordinates for each group
    coords_a = all_coords[indices_a]
    coords_b = all_coords[indices_b]

    # Compute centers of mass
    com_a = compute_com(coords_a)
    com_b = compute_com(coords_b)

    print(f"\nCOM of A: ({com_a[0]:.6f}, {com_a[1]:.6f}, {com_a[2]:.6f})")
    print(f"COM of B: ({com_b[0]:.6f}, {com_b[1]:.6f}, {com_b[2]:.6f})")

    # Distance between COMs
    distance = np.linalg.norm(com_a - com_b)
    print(f"\nDistance between COMs: {distance:.6f}")

if __name__ == "__main__":
    main()

