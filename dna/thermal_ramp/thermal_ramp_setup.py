#!/usr/bin/env python3

import argparse
import os
import glob
import shutil
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description="Setup thermal ramp MD inputs for oxDNA.")
    parser.add_argument("-s", type=int, required=True, help="Starting temperature (C)")
    parser.add_argument("-e", type=int, required=True, help="Ending temperature (C)")
    parser.add_argument("-t", type=int, required=True, help="Number of ramp steps")
    parser.add_argument("-input", required=True, help="Reference mdinput file")
    return parser.parse_args()


def find_relax_dir():
    relax_dirs = glob.glob("*_relax")

    if len(relax_dirs) == 0:
        sys.exit("ERROR: No *_relax directory found.")

    if len(relax_dirs) > 1:
        sys.exit("ERROR: Multiple *_relax directories found. Keep only one.")

    return Path(relax_dirs[0])


def find_top_file(relax_dir):
    tops = list(relax_dir.glob("*.top"))

    if len(tops) == 0:
        sys.exit("ERROR: No .top file found in relax directory.")

    if len(tops) > 1:
        sys.exit("ERROR: Multiple .top files found.")

    return tops[0]


def load_template(path):
    if not os.path.exists(path):
        sys.exit(f"ERROR: input template {path} not found.")

    with open(path) as f:
        return f.readlines()


def linear_ramp(start, end, steps):

    if steps == 1:
        return [start]

    delta = (end - start) / (steps - 1)

    return [round(start + i * delta) for i in range(steps)]


def replace_param(lines, key, value):

    new = []

    for line in lines:

        if line.strip().startswith(key):
            new.append(f"{key} = {value}\n")

        else:
            new.append(line)

    return new


def write_input(path, lines):

    with open(path, "w") as f:
        f.writelines(lines)


def main():

    args = parse_args()

    start = args.s
    end = args.e
    steps = args.t
    template = args.input

    relax_dir = find_relax_dir()
    prefix = relax_dir.name.replace("_relax", "")

    ramp_dir = Path(f"{prefix}_thermal_ramp")

    if ramp_dir.exists():
        sys.exit(f"ERROR: {ramp_dir} already exists.")

    ramp_dir.mkdir()

    top_file = find_top_file(relax_dir)
    conf_file = relax_dir / "md_relax_last_conf.dat"

    if not conf_file.exists():
        sys.exit("ERROR: md_relax_last_conf.dat not found.")

    print(f"Using relax directory: {relax_dir}")
    print(f"Thermal ramp directory: {ramp_dir}")
    print(f"Topology: {top_file}")
    print(f"Initial configuration: {conf_file}")

    temps = linear_ramp(start, end, steps)

    template_lines = load_template(template)

    prev_conf = conf_file.name

    run_commands = []
    summary = []

    for i, T in enumerate(temps, start=1):

        step_name = f"step{i}_{T}c"

        input_name = f"inputMD_{step_name}"   # <-- NO SUFFIX

        conf_out = f"{step_name}_last_conf.dat"
        traj_out = f"{step_name}_trajectory.dat"
        energy_out = f"{step_name}_energy.dat"

        lines = template_lines.copy()

        lines = replace_param(lines, "topology", top_file.name)
        lines = replace_param(lines, "conf_file", prev_conf)
        lines = replace_param(lines, "T", f"{T}C")
        lines = replace_param(lines, "lastconf_file", conf_out)
        lines = replace_param(lines, "trajectory_file", traj_out)
        lines = replace_param(lines, "energy_file", energy_out)

        if i == 1:
            lines = replace_param(lines, "restart_step_counter", "1")
        else:
            lines = replace_param(lines, "restart_step_counter", "0")

        input_path = ramp_dir / input_name

        write_input(input_path, lines)

        run_commands.append(f"/home/dfan/oxDNA/build/bin/oxDNA {input_name}")
        run_commands.append("sleep 5")

        summary.append((i, T, prev_conf, conf_out))

        prev_conf = conf_out

    shutil.copy(top_file, ramp_dir)
    shutil.copy(conf_file, ramp_dir / conf_file.name)

    run_script = ramp_dir / "run_thermal_ramp.sh"

    with open(run_script, "w") as f:

        f.write("#!/bin/bash\n\n")

        for cmd in run_commands:
            f.write(cmd + "\n")

    os.chmod(run_script, 0o755)

    print("\n================ SUMMARY ================\n")

    print("Step | Temperature | Input Conf -> Output Conf")
    print("----------------------------------------------")

    for s in summary:
        print(f"{s[0]:>4} | {s[1]:>10}C | {s[2]} -> {s[3]}")

    print("\nRun script:")
    print(run_script)

    print("\n========== THERMAL RAMP SETTINGS ==========")
    print(f"Start Temperature : {start}C")
    print(f"End Temperature   : {end}C")
    print(f"Number of steps   : {steps}")

    print("\nRamp sequence:")

    for i, T in enumerate(temps, start=1):
        print(f"Step {i} → {T}C")

    print("\n===========================================")


if __name__ == "__main__":
    main()
