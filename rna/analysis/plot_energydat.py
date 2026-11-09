#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import os
import sys

def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot oxDNA potential energy vs timestep."
    )
    parser.add_argument(
        "energy_file",
        help="energy.dat file to plot"
    )
    return parser.parse_args()


def read_energy_file(filename):

    timesteps = []
    energies = []

    with open(filename) as f:
        for line in f:

            if not line.strip():
                continue

            parts = line.split()

            if len(parts) < 2:
                continue

            try:
                # Remove decimal point completely and convert to int
                timestep_str = parts[0].replace('.', '')  # '19990.0000' -> '19990000'
                timestep = int(timestep_str)

                potential_energy = float(parts[1])

                timesteps.append(timestep)
                energies.append(potential_energy)

            except ValueError:
                continue

    if len(timesteps) == 0:
        print("ERROR: No valid energy data found.")
        sys.exit(1)

    return timesteps, energies


def main():

    args = parse_args()

    if not os.path.exists(args.energy_file):
        print(f"ERROR: File not found: {args.energy_file}")
        sys.exit(1)

    print(f"Reading energy file: {args.energy_file}")

    timesteps, energies = read_energy_file(args.energy_file)

    plt.figure(figsize=(8, 5))

    plt.plot(timesteps, energies)

    plt.xlabel("Timestep")
    plt.ylabel("Potential Energy")
    plt.title("Potential Energy vs Timestep")

    plt.tight_layout()

    output_name = "energy_plot.png"

    plt.savefig(output_name, dpi=300)

    print(f"Plot saved to: {output_name}")


if __name__ == "__main__":
    main()
