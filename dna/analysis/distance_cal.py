#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np

def calculate_distance(p1, p2):
    return np.linalg.norm(np.array(p1) - np.array(p2))

def main():

    if len(sys.argv) != 4:
        print("Usage: python distance_calculator.py <trj.dat> <nt1> <nt2>")
        sys.exit(1)

    traj_file = sys.argv[1]
    nt1 = int(sys.argv[2])
    nt2 = int(sys.argv[3])

    timestep = None
    particle_index = 0
    coords = {}
    results = []

    with open(traj_file, "r") as f:
        for line in f:
            line = line.strip()

            if line.startswith("t ="):
                if coords:
                    if nt1 in coords and nt2 in coords:
                        dist = calculate_distance(coords[nt1], coords[nt2])
                        results.append((timestep, dist))

                timestep = int(line.split("=")[1])
                coords = {}
                particle_index = 0
                continue

            if line.startswith("b =") or line.startswith("E =") or line == "":
                continue

            parts = line.split()

            if len(parts) >= 3:
                particle_index += 1
                x, y, z = map(float, parts[:3])
                coords[particle_index] = (x, y, z)

        # handle last timestep
        if coords and nt1 in coords and nt2 in coords:
            dist = calculate_distance(coords[nt1], coords[nt2])
            results.append((timestep, dist))

    df = pd.DataFrame(results, columns=["timestep", "distance"])

    output_name = f"{nt1}-{nt2}_distance.csv"
    df.to_csv(output_name, index=False)

    print(f"Saved distance data to {output_name}")


if __name__ == "__main__":
    main()
