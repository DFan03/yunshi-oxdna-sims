#!/usr/bin/env python3

import os
import glob
import shutil
import sys


def parse_input_parameters(filepath):
    params = {
        "T": None,
        "use_average_seq": None,
        "dt": None,
        "mismatch_repulsion": None,
        "interaction_type": None
    }

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()

            if "=" not in line:
                continue

            key, value = [x.strip() for x in line.split("=", 1)]

            if key in params:
                params[key] = value

    # normalize booleans
    for k in ["use_average_seq", "mismatch_repulsion"]:
        v = params[k]
        if v is None:
            continue

        v = v.lower()
        if v in ["1", "true"]:
            params[k] = "True"
        elif v in ["0", "false"]:
            params[k] = "False"

    return params


def main():

    cwd = os.getcwd()
    prefix = os.path.basename(cwd)

    # --------------------------------------------------
    # Detect files
    # --------------------------------------------------
    top_files = glob.glob("*.top")
    oxdna_files = glob.glob("*.oxdna")
    dat_files = glob.glob("*.dat")

    if len(top_files) != 1:
        print(f"ERROR: Expected exactly 1 .top file, found {len(top_files)}")
        sys.exit(1)

    if oxdna_files and dat_files:
        print("ERROR: Both .oxdna and .dat files found. Only one allowed.")
        sys.exit(1)

    if not oxdna_files and not dat_files:
        print("ERROR: No structure file (.oxdna or .dat) found.")
        sys.exit(1)

    structure_files = oxdna_files + dat_files

    if len(structure_files) != 1:
        print(f"ERROR: Expected exactly 1 structure file, found {len(structure_files)}")
        sys.exit(1)

    top_file = top_files[0]
    structure_file = structure_files[0]

    # --------------------------------------------------
    # Define workflow directories
    # --------------------------------------------------
    relax_dir = f"{prefix}_relax"
    md_dir = f"{prefix}_md"
    analysis_dir = f"{prefix}_analysis"

    workflow_dirs = [relax_dir, md_dir, analysis_dir]

    # --------------------------------------------------
    # Prevent overwriting
    # --------------------------------------------------
    for d in workflow_dirs:
        if os.path.exists(d):
            print(f"ERROR: Directory '{d}' already exists.")
            print("Refusing to overwrite existing workflow directories.")
            sys.exit(1)

    # --------------------------------------------------
    # Create directories
    # --------------------------------------------------
    for d in workflow_dirs:
        os.makedirs(d)

    # --------------------------------------------------
    # Copy core files
    # --------------------------------------------------
    for d in workflow_dirs:
        shutil.copy(top_file, d)

    shutil.copy(structure_file, relax_dir)

    # --------------------------------------------------
    # Copy utilities
    # --------------------------------------------------
    utils_dir = "/home/dfan/rna_utils"

    shutil.copy(os.path.join(utils_dir, "input_pre_relax"), relax_dir)
    shutil.copy(os.path.join(utils_dir, "input_cuda_relax"), relax_dir)
    shutil.copy(os.path.join(utils_dir, "plot_energy.ipynb"), relax_dir)


    # --------------------------------------------------
    # Replace placeholders
    # --------------------------------------------------
    def replace_placeholders(filepath):
        with open(filepath, "r") as f:
            lines = f.readlines()

        new_lines = []
        for line in lines:
            line = line.replace("replace_me.top", top_file)
            line = line.replace("replace_me.dat", structure_file)
            new_lines.append(line)

        with open(filepath, "w") as f:
            f.writelines(new_lines)

    pre_relax_path = os.path.join(relax_dir, "input_pre_relax")
    cuda_relax_path = os.path.join(relax_dir, "input_cuda_relax")

    replace_placeholders(pre_relax_path)
    replace_placeholders(cuda_relax_path)

    # --------------------------------------------------
    # Create cmds_relax
    # --------------------------------------------------
    cmds_path = os.path.join(relax_dir, "cmds_relax")

    with open(cmds_path, "w") as f:
        f.write("/home/dfan/oxDNA/build/bin/oxDNA input_pre_relax\n")
        f.write("sleep 5\n")
        f.write("/home/dfan/oxDNA/build/bin/oxDNA input_cuda_relax\n")

    # --------------------------------------------------
    # Parse configuration parameters
    # --------------------------------------------------
    mc_params = parse_input_parameters(pre_relax_path)
    md_params = parse_input_parameters(cuda_relax_path)

    # --------------------------------------------------
    # Print summary
    # --------------------------------------------------
    print("\nWorkflow setup complete\n")

    print("Detected files")
    print("---------------------------")
    print(f"Topology   : {top_file}")
    print(f"Structure  : {structure_file}")

    print("\nDirectories created")
    print("---------------------------")
    for d in workflow_dirs:
        print(d)

    print("\nFiles placed")
    print("---------------------------")

    print(f"\n{relax_dir}")
    print(f"  {top_file}")
    print(f"  {structure_file}")
    print("  input_pre_relax")
    print("  input_cuda_relax")
    print("  plot_energy.ipynb")
    print("  cmds_relax")

    print(f"\n{md_dir}")
    print(f"  {top_file}")

    print(f"\n{analysis_dir}")
    print(f"  {top_file}")
    print("  prepare_oat_workflow.py")
    print("  calculate_bond_instability.py")

    print("\nSimulation configuration")
    print("---------------------------")

    print("\nMC relax configuration (input_pre_relax)")
    print(f"  Temperature (T)          : {mc_params['T']}")
    print(f"  use_average_seq          : {mc_params['use_average_seq']}")
    print(f"  dt                       : {mc_params['dt']}")
    print(f"  mismatch_repulsion       : {mc_params['mismatch_repulsion']}")
    print(f"  interaction_type         : {mc_params['interaction_type']}")

    print("\nMD relax configuration (input_cuda_relax)")
    print(f"  Temperature (T)          : {md_params['T']}")
    print(f"  use_average_seq          : {md_params['use_average_seq']}")
    print(f"  dt                       : {md_params['dt']}")
    print(f"  mismatch_repulsion       : {md_params['mismatch_repulsion']}")
    print(f"  interaction_type         : {md_params['interaction_type']}")

    print("\nDone.")


if __name__ == "__main__":
    main()
