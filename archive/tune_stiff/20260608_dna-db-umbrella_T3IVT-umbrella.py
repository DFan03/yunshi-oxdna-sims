#!/usr/bin/env python3
"""
tune_stiff_prod.py

Usage:
    python tune_stiff_prod.py NEW_STIFF NEW_STEPS NEW_PRINT_CONF NEW_PRINT_ENERGY [OXDNA_BIN]

Run this script from inside the umbrella directory (where window_00, window_01, ... live).

What it does:
  - For each window_XX/prod/forces.txt:
        * replaces all 'stiff = ...' lines with 'stiff = NEW_STIFF'
  - For each window_XX/prod/input.txt:
        * sets  'steps = NEW_STEPS'
        * sets  'print_conf_interval = NEW_PRINT_CONF'
        * sets  'print_energy_every  = NEW_PRINT_ENERGY'
  - Deletes old production outputs in each prod/:
        * trajectory.dat, last_conf.dat, energy.dat, run.log (if present)
  - Writes a new run_production_only.sh to re-run production for all windows
        * uses OXDNA_BIN if provided, else defaults to 'oxDNA'
"""

import os
import sys
import glob

# -----------------------------
# Helpers
# -----------------------------

def update_stiff_in_forces(forces_path, new_stiff):
    """Replace all 'stiff = ...' lines in the given forces.txt with new_stiff."""
    print(f"[update_stiff_in_forces] Updating {forces_path}")
    if not os.path.isfile(forces_path):
        raise FileNotFoundError(f"forces.txt not found: {forces_path}")

    with open(forces_path) as f:
        lines = f.readlines()

    new_lines = []
    replacements = 0
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("stiff"):
            # Preserve indentation
            indent = line[:len(line) - len(line.lstrip())]
            new_line = f"{indent}stiff = {new_stiff}\n"
            new_lines.append(new_line)
            replacements += 1
        else:
            new_lines.append(line)

    if replacements == 0:
        print(f"  [warn] No 'stiff' lines found in {forces_path} (nothing changed).")
    else:
        print(f"  Replaced {replacements} 'stiff' lines with stiff = {new_stiff}")

    with open(forces_path, 'w') as f:
        f.writelines(new_lines)


def update_prod_input(input_path, new_steps, new_print_conf, new_print_energy):
    """
    Update 'steps', 'print_conf_interval', and 'print_energy_every'
    in a production input.txt, using user-specified values.
    """
    print(f"[update_prod_input] Updating {input_path}")
    if not os.path.isfile(input_path):
        raise FileNotFoundError(f"input.txt not found: {input_path}")

    with open(input_path) as f:
        lines = f.readlines()

    old_steps = None
    old_print_conf = None
    old_print_energy = None

    # Read old values for logging
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("steps"):
            parts = stripped.replace(" ", "").split("=")
            if len(parts) == 2:
                try:
                    old_steps = int(parts[1])
                except ValueError:
                    pass
        elif stripped.startswith("print_conf_interval"):
            parts = stripped.replace(" ", "").split("=")
            if len(parts) == 2:
                try:
                    old_print_conf = int(parts[1])
                except ValueError:
                    pass
        elif stripped.startswith("print_energy_every"):
            parts = stripped.replace(" ", "").split("=")
            if len(parts) == 2:
                try:
                    old_print_energy = int(parts[1])
                except ValueError:
                    pass

    new_lines = []
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("steps"):
            indent = line[:len(line) - len(line.lstrip())]
            new_lines.append(f"{indent}steps = {new_steps}\n")
        elif stripped.startswith("print_conf_interval"):
            indent = line[:len(line) - len(line.lstrip())]
            new_lines.append(f"{indent}print_conf_interval = {new_print_conf}\n")
        elif stripped.startswith("print_energy_every"):
            indent = line[:len(line) - len(line.lstrip())]
            new_lines.append(f"{indent}print_energy_every = {new_print_energy}\n")
        else:
            new_lines.append(line)

    with open(input_path, 'w') as f:
        f.writelines(new_lines)

    print("  Updated input.txt:")
    if old_steps is not None:
        print(f"    steps              : {old_steps} -> {new_steps}")
    else:
        print(f"    steps              : (set to {new_steps})")
    if old_print_conf is not None:
        print(f"    print_conf_interval: {old_print_conf} -> {new_print_conf}")
    else:
        print(f"    print_conf_interval: (set to {new_print_conf})")
    if old_print_energy is not None:
        print(f"    print_energy_every : {old_print_energy} -> {new_print_energy}")
    else:
        print(f"    print_energy_every : (set to {new_print_energy})")


def remove_old_production_files(prod_dir):
    """
    Remove old production outputs from prod_dir.
    Does NOT touch equil directories.
    """
    print(f"[remove_old_production_files] Cleaning {prod_dir}")
    targets = ["trajectory.dat", "last_conf.dat", "energy.dat", "run.log"]
    for name in targets:
        path = os.path.join(prod_dir, name)
        if os.path.isfile(path):
            os.remove(path)
            print(f"  Removed {path}")
        else:
            print(f"  (no {name} to remove)")


def find_window_dirs(pattern="window_*"):
    """Return sorted list of window directories (window_00, window_01, ...)."""
    dirs = sorted(d for d in glob.glob(pattern) if os.path.isdir(d))
    return dirs


def write_run_script(window_dirs, oxdna_bin, script_name="run_production_only.sh"):
    """
    Create a shell script that runs production in all windows sequentially.
    """
    print(f"[write_run_script] Writing {script_name}")
    lines = []
    lines.append("#!/bin/bash")
    lines.append("# Auto-generated by tune_stiff_prod.py")
    lines.append("set -e")
    lines.append("")

    for w in window_dirs:
        prod_dir = os.path.join(w, "prod")
        lines.append(f"cd {prod_dir} && {oxdna_bin} input.txt > run.log 2>&1 && cd ..")

    script_path = os.path.join(os.getcwd(), script_name)
    with open(script_path, 'w') as f:
        f.write("\n".join(lines) + "\n")

    os.chmod(script_path, 0o755)
    print(f"[write_run_script] Created {script_name}")


# -----------------------------
# Main
# -----------------------------

def main():
    if len(sys.argv) < 5:
        print("Usage: python tune_stiff_prod.py NEW_STIFF NEW_STEPS NEW_PRINT_CONF NEW_PRINT_ENERGY [OXDNA_BIN]")
        sys.exit(1)

    try:
        new_stiff = float(sys.argv[1])
    except ValueError:
        print(f"ERROR: NEW_STIFF must be a number, got {sys.argv[1]!r}")
        sys.exit(1)

    try:
        new_steps = int(sys.argv[2])
        if new_steps <= 0:
            raise ValueError
    except ValueError:
        print(f"ERROR: NEW_STEPS must be a positive integer, got {sys.argv[2]!r}")
        sys.exit(1)

    try:
        new_print_conf = int(sys.argv[3])
        if new_print_conf <= 0:
            raise ValueError
    except ValueError:
        print(f"ERROR: NEW_PRINT_CONF must be a positive integer, got {sys.argv[3]!r}")
        sys.exit(1)

    try:
        new_print_energy = int(sys.argv[4])
        if new_print_energy <= 0:
            raise ValueError
    except ValueError:
        print(f"ERROR: NEW_PRINT_ENERGY must be a positive integer, got {sys.argv[4]!r}")
        sys.exit(1)

    if len(sys.argv) >= 6:
        oxdna_bin = sys.argv[5]
    else:
        oxdna_bin = "oxDNA"  # default; edit if needed

    print("=== tune_stiff_prod.py ===")
    print(f"Current directory        : {os.getcwd()}")
    print(f"New stiffness            : {new_stiff}")
    print(f"New production steps     : {new_steps}")
    print(f"New print_conf_interval  : {new_print_conf}")
    print(f"New print_energy_every   : {new_print_energy}")
    print(f"oxDNA binary             : {oxdna_bin}")
    print("")

    # Find all windows
    window_dirs = find_window_dirs("window_*")
    if not window_dirs:
        print("ERROR: No window_* directories found here. Run from the umbrella directory.")
        sys.exit(1)

    print(f"Found {len(window_dirs)} window(s):")
    for w in window_dirs:
        print(f"  - {w}")
    print("")

    # Update forces/input and clean prod dirs
    for w in window_dirs:
        print(f"\n[Window] {w}")
        prod_dir = os.path.join(w, "prod")
        forces_path = os.path.join(prod_dir, "forces.txt")
        input_path  = os.path.join(prod_dir, "input.txt")

        if not os.path.isdir(prod_dir):
            print(f"  WARNING: {prod_dir} does not exist, skipping.")
            continue

        # 1) Update stiffness in forces.txt
        update_stiff_in_forces(forces_path, new_stiff)

        # 2) Update steps and intervals in input.txt
        update_prod_input(input_path, new_steps, new_print_conf, new_print_energy)

        # 3) Remove old production outputs
        remove_old_production_files(prod_dir)

    # 4) Create a new run script for production only
    write_run_script(window_dirs, oxdna_bin, script_name="run_production_only.sh")

    print("\nAll done.")
    print("Run new production with:")
    print("  bash run_production_only.sh")


if __name__ == "__main__":
    main()
