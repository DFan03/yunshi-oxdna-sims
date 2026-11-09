#!/usr/bin/env python3

import sys
import pandas as pd
import matplotlib.pyplot as plt
import os


def parse_energy_file(energy_file):
    timesteps = []
    energies = []

    with open(energy_file, "r") as f:
        for line in f:
            parts = line.split()

            if len(parts) < 2:
                continue

            # convert timestep
            # remove decimal part and append 000
            t_raw = parts[0]
            t = int(t_raw.split(".")[0]) * 1000

            energy = float(parts[3])

            timesteps.append(t)
            energies.append(energy)

    df = pd.DataFrame({
        "timestep": timesteps,
        "energy": energies
    })

    return df


def main():

    if len(sys.argv) < 2:
        print("Usage:")
        print("dis_energy_plotter <distance.csv> <energy.dat> <dis_only=false>")
        sys.exit(1)

    distance_file = sys.argv[1]

    energy_file = None
    dis_only = False

    if len(sys.argv) >= 3:
        energy_file = sys.argv[2]

    if len(sys.argv) == 4:
        dis_only = sys.argv[3].lower() == "true"

    # read distance file
    dist_df = pd.read_csv(distance_file)

    distance_name = os.path.basename(distance_file).replace(".csv", "")

    if dis_only:

        plt.figure(figsize=(30,5))

        plt.plot(
            dist_df["timestep"],
            dist_df["distance"],
            color="tab:blue",
            alpha=0.5,
            linewidth=2,
            label=distance_name
        )

        plt.xlabel("t")
        plt.ylabel("distance")
        plt.legend()

        plt.tight_layout()
        plt.show()

        return

    if energy_file is None:
        print("Energy file required unless dis_only=true")
        sys.exit(1)

    energy_df = parse_energy_file(energy_file)

    # merge datasets by timestep
    merged = pd.merge(dist_df, energy_df, on="timestep", how="inner")

    # create elongated plot
    fig, ax1 = plt.subplots(figsize=(14,5))

    # energy axis (left)
    ax1.set_xlabel("t")
    ax1.set_ylabel("energy")

    energy_line = ax1.plot(
        merged["timestep"],
        merged["energy"],
        color="tab:red",
        alpha=0.5,
        linewidth=1,
        label="potential energy"
    )

    # distance axis (right)
    ax2 = ax1.twinx()
    ax2.set_ylabel("distance")

    dist_line = ax2.plot(
        merged["timestep"],
        merged["distance"],
        color="tab:blue",
        alpha=0.5,
        linewidth=1,
        label=distance_name
    )

    # combine legends
    lines = energy_line + dist_line
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, loc="best")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
