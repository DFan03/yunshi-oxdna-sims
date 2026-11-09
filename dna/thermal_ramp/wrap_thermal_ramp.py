#!/usr/bin/env python3

import os
import glob
import shutil
import subprocess
import re
import sys


def natural_sort_key(s):
    """
    Sort filenames like step1, step2, step10 correctly.
    """
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(r'(\d+)', s)]


def move_files(pattern, destination):

    files = glob.glob(pattern)

    if not files:
        print(f"[WARNING] No files found for pattern: {pattern}")
        return []

    os.makedirs(destination, exist_ok=True)

    moved = []

    for f in files:
        dest = os.path.join(destination, os.path.basename(f))
        shutil.move(f, dest)
        moved.append(dest)

    return moved


def concatenate(files, output):

    if not files:
        print(f"[WARNING] No files to concatenate for {output}")
        return

    with open(output, "w") as outfile:
        for f in files:
            with open(f) as infile:
                outfile.write(infile.read())

    print(f"[OK] Created {output}")


def main():

    print("\n=== oxDNA Output Organizer ===\n")

    energy_files = glob.glob("*_energy.dat")
    traj_files = glob.glob("*_trajectory.dat")
    conf_files = glob.glob("*_last_conf.dat")
    input_files = glob.glob("inputMD_step*")

    if not (energy_files or traj_files or conf_files or input_files):
        print("ERROR: No oxDNA related files detected.")
        sys.exit(1)

    print("Creating directories and moving files...\n")

    energy_moved = move_files("*_energy.dat", "energy")
    traj_moved = move_files("*_trajectory.dat", "traj")
    conf_moved = move_files("*_last_conf.dat", "last_conf")
    input_moved = move_files("inputMD_step*", "inputs")

    # ensure correct numerical order
    energy_moved.sort(key=natural_sort_key)
    traj_moved.sort(key=natural_sort_key)

    print("\nConcatenating energy files...")
    concatenate(energy_moved, "energy_total.dat")

    print("\nConcatenating trajectory files...")
    concatenate(traj_moved, "traj_total.dat")

    print("\nRunning oat trajectory alignment...")

    try:
        subprocess.run(
            ["oat", "align", "traj_total.dat", "aligned_trj.dat"],
            check=True
        )
        print("[OK] aligned_trj.dat generated")

    except Exception as e:
        print("[WARNING] oat alignment failed.")
        print(e)

    print("\n===== SUMMARY =====")

    print(f"Energy files moved : {len(energy_moved)}")
    print(f"Trajectory files   : {len(traj_moved)}")
    print(f"Last conf files    : {len(conf_moved)}")
    print(f"Input files moved  : {len(input_moved)}")

    print("\nDirectories created:")
    print("energy/")
    print("traj/")
    print("last_conf/")
    print("inputs/")

    print("\nGenerated files:")
    print("energy_total.dat")
    print("traj_total.dat")
    print("aligned_trj.dat")

    print("\n====================\n")


if __name__ == "__main__":
    main()
