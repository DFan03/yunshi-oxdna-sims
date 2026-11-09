#!/usr/bin/env python3
"""
setup_kl_ffs.py

Usage:
    python setup_kl_ffs.py X

Run this in a directory that contains exactly ONE *.top and ONE *.dat file.
It will:

- Copy /home/dfan/rna_utils/fss/op.txt here and read how many bond pairs there are.
- Create:
    FLUX,
    and shooting dirs: (N-1) down to 1, like 54_shoot, 43_shoot, ..., 10_shoot.
- Copy templates and configure:
    - FLUX/ffs_flux.py
    - each */input
    - each */ffs_shoot.py
- Create a cmd_ffs bash script in the current directory that runs all stages.
"""

import os
import sys
import glob
import shutil
import re

# ---------------------------------------------------------------------
# User-adjustable parameters
# ---------------------------------------------------------------------

TEMPLATE_DIR = "/home/dfan/rna_utils/fss"

# Default number of MD steps
FLUX_STEPS = "100"   # steps in FLUX/input
SHOOT_STEPS = "2e10"  # steps in */input for all shooting stages

NCPUS = 8             # number of CPUs used in cmd_ffs


# ---------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------


def find_single(pattern: str, description: str) -> str:
    matches = glob.glob(pattern)
    if len(matches) == 0:
        raise SystemExit(f"Error: no {description} matching {pattern} found in {os.getcwd()}")
    if len(matches) > 1:
        raise SystemExit(
            f"Error: multiple {description} matching {pattern} found:\n  " + "\n  ".join(matches)
        )
    return matches[0]


def copy_if_not_exists(src: str, dst: str) -> None:
    if not os.path.exists(src):
        raise SystemExit(f"Error: template file not found: {src}")
    shutil.copy(src, dst)


def parse_pair_count_from_op(op_path: str) -> int:
    n_pairs = 0
    with open(op_path, "r") as f:
        for line in f:
            line_stripped = line.strip()
            if line_stripped.startswith("pair"):
                n_pairs += 1
    if n_pairs == 0:
        raise SystemExit(f"Error: no 'pairX =' lines found in {op_path}")
    return n_pairs


def replace_in_input(
    input_path: str,
    top_name: str,
    conf_name: str | None = None,
    steps_val: str | None = None,
) -> None:
    """
    Replace topology=, optionally conf_file=, and optionally steps= in an oxDNA input file.
    """
    new_lines = []
    with open(input_path, "r") as f:
        for line in f:
            stripped = line.lstrip()
            indent = line[: len(line) - len(stripped)]
            if stripped.startswith("topology"):
                new_lines.append(f"{indent}topology = {top_name}\n")
            elif conf_name is not None and stripped.startswith("conf_file"):
                new_lines.append(f"{indent}conf_file = {conf_name}\n")
            elif steps_val is not None and stripped.startswith("steps"):
                new_lines.append(f"{indent}steps = {steps_val}\n")
            else:
                new_lines.append(line)

    with open(input_path, "w") as f:
        f.writelines(new_lines)


def configure_ffs_flux(flux_dir: str, dat_name: str, template_flux: str) -> None:
    dst = os.path.join(flux_dir, "ffs_flux.py")
    copy_if_not_exists(template_flux, dst)

    with open(dst, "r") as f:
        txt = f.read()

    txt_new, n_sub = re.subn(
        r"starting_conf\s*=\s*['\"].*?['\"]",
        f"starting_conf = '{dat_name}'",
        txt,
    )
    if n_sub == 0:
        raise SystemExit(
            f"Error: could not find a line like \"starting_conf = '...'\") in {dst}"
        )

    with open(dst, "w") as f:
        f.write(txt_new)


def configure_ffs_shoot(
    shoot_dir: str,
    top_name: str,
    upstream_dir: str,
    high_val: int,
    low_val: int,
    template_shoot: str,
) -> None:
    dst = os.path.join(shoot_dir, "ffs_shoot.py")
    copy_if_not_exists(template_shoot, dst)

    with open(dst, "r") as f:
        lines = f.readlines()

    new_lines = []
    for line in lines:
        stripped = line.lstrip()
        indent = line[: len(line) - len(stripped)]

        if stripped.startswith("lambda_f_value"):
            new_lines.append(f"{indent}lambda_f_value = {high_val}.\n")
        elif stripped.startswith("lambda_n_value"):
            new_lines.append(f"{indent}lambda_n_value = {low_val}.\n")
        elif stripped.startswith("lambda_m_value"):
            new_lines.append(f"{indent}lambda_m_value = {low_val}\n")
        elif stripped.startswith("starting_conf_pattern"):
            pattern = f"../{upstream_dir}/success*"
            new_lines.append(f"{indent}starting_conf_pattern = '{pattern}'\n")
        elif stripped.startswith("topology"):
            new_lines.append(f"{indent}topology = '{top_name}'\n")
        else:
            new_lines.append(line)

    with open(dst, "w") as f:
        f.writelines(new_lines)


def parse_input_params_for_report(input_path: str) -> dict:
    keys = ["T", "use_average_seq", "interaction_type", "mismatch_repulsion", "steps"]
    values = {k: None for k in keys}

    with open(input_path, "r") as f:
        for line in f:
            stripped = line.strip()
            if "=" not in stripped:
                continue
            left, right = stripped.split("=", 1)
            key = left.strip()
            if key in values:
                values[key] = right.strip()

    return values


def make_cmd_script(
    script_path: str,
    target_successes: int,
    flux_dir: str,
    shoot_dirs: list[str],
) -> None:
    lines = [
        "#!/usr/bin/env bash",
        "set -e",
        "",
        f"X={target_successes}",
        "",
        'echo "Running FLUX stage..."',
        f"( cd {flux_dir} && python ffs_flux.py -c {NCPUS} -n \"$X\" -v )",
        "sleep 2",
        "",
    ]

    for d in shoot_dirs:
        lines.append(f'echo "Running shooting stage in {d}..."')
        lines.append(f"( cd {d} && python ffs_shoot.py -c {NCPUS} -n \"$X\" -v )")
        lines.append("sleep 2")
        lines.append("")

    with open(script_path, "w") as f:
        f.write("\n".join(lines) + "\n")

    os.chmod(script_path, 0o755)


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------


def main():
    if len(sys.argv) != 2:
        print("Usage: python setup_kl_ffs.py X", file=sys.stderr)
        sys.exit(1)

    try:
        target_successes = int(sys.argv[1])
    except ValueError:
        raise SystemExit("Error: X must be an integer (desired number of successes)")

    cwd = os.getcwd()

    # 1. Find top and dat files in current directory
    top_file = find_single("*.top", "*.top file")
    dat_file = find_single("*.dat", "*.dat file")

    # 2. Template paths
    op_src = os.path.join(TEMPLATE_DIR, "op.txt")
    input_src = os.path.join(TEMPLATE_DIR, "input")
    flux_src = os.path.join(TEMPLATE_DIR, "ffs_flux.py")
    shoot_src = os.path.join(TEMPLATE_DIR, "ffs_shoot.py")

    # Copy op.txt to cwd
    op_dst = os.path.join(cwd, "op.txt")
    copy_if_not_exists(op_src, op_dst)

    # 3. Parse number of pairs from op.txt
    n_pairs = parse_pair_count_from_op(op_dst)
    if n_pairs < 1:
        raise SystemExit("Error: parsed number of pairs < 1, this should not happen.")

    # 4. Build directory list
    flux_dir = "FLUX"
    shoot_dirs: list[str] = []
    for i in range(n_pairs - 1, 0, -1):
        shoot_dirs.append(f"{i}{i-1}_shoot")

    all_dirs = [flux_dir] + shoot_dirs

    # 5. Create directories and populate basic files
    for d in all_dirs:
        if os.path.exists(d):
            raise SystemExit(f"Error: directory {d} already exists. Aborting to avoid overwrite.")
        os.mkdir(d)

        shutil.copy(top_file, os.path.join(d, os.path.basename(top_file)))
        shutil.copy(op_dst, os.path.join(d, "op.txt"))
        copy_if_not_exists(input_src, os.path.join(d, "input"))

    # .dat only into FLUX
    shutil.copy(dat_file, os.path.join(flux_dir, os.path.basename(dat_file)))

    # 6. Configure FLUX
    configure_ffs_flux(flux_dir, os.path.basename(dat_file), flux_src)

    flux_input_path = os.path.join(flux_dir, "input")
    replace_in_input(
        flux_input_path,
        top_name=os.path.basename(top_file),
        conf_name=os.path.basename(dat_file),
        steps_val=FLUX_STEPS,
    )

    # 7. Configure shooting directories
    top_name_only = os.path.basename(top_file)

    for idx, d in enumerate(shoot_dirs):
        prefix = d.split("_")[0]
        if len(prefix) != 2 or not prefix.isdigit():
            raise SystemExit(
                f"Error: unexpected shoot directory name '{d}'. "
                "Expected a two-digit prefix like '54_shoot'."
            )
        high_val = int(prefix[0])
        low_val = int(prefix[1])
        upstream = flux_dir if idx == 0 else shoot_dirs[idx - 1]

        configure_ffs_shoot(
            shoot_dir=d,
            top_name=top_name_only,
            upstream_dir=upstream,
            high_val=high_val,
            low_val=low_val,
            template_shoot=shoot_src,
        )

        shoot_input_path = os.path.join(d, "input")
        replace_in_input(
            shoot_input_path,
            top_name=top_name_only,
            conf_name=None,
            steps_val=SHOOT_STEPS,
        )

    # 8. Create cmd_ffs script in master directory
    cmd_script_path = os.path.join(cwd, "cmd_ffs")
    make_cmd_script(
        script_path=cmd_script_path,
        target_successes=target_successes,
        flux_dir=flux_dir,
        shoot_dirs=shoot_dirs,
    )

    # 9. Read input parameters for report (from FLUX/input)
    params = parse_input_params_for_report(flux_input_path)

    # 10. Print summary
    print("\n=== FFS setup complete ===\n")
    print(f"Working directory: {cwd}")
    print(f"Top file:          {top_file}")
    print(f"Conf file (.dat):  {dat_file}")
    print(f"Number of pairs:   {n_pairs}")
    print(f"Target successes X (for cmd_ffs): {target_successes}")
    print()
    print("Created directories (in order):")
    for d in all_dirs:
        print(f"  - {d}")
    print()
    print("Run all stages with:")
    print("  bash cmd_ffs")
    print()
    print("Key input parameters (from FLUX/input):")
    for k in ["T", "use_average_seq", "interaction_type", "mismatch_repulsion", "steps"]:
        print(f"  {k} = {params.get(k)}")
    print()
    print("Summary table:")
    print("-----------------------------------------------")
    print(f"{'Stage':<15} {'Upstream':<15} {'lambda_f -> lambda_m':<15}")
    print("-----------------------------------------------")
    print(f"{'FLUX':<15} {'-':<15} {f'{n_pairs} -> {n_pairs-1}':<15}")
    for idx, d in enumerate(shoot_dirs):
        prefix = d.split("_")[0]
        high_val = int(prefix[0])
        low_val = int(prefix[1])
        upstream = "FLUX" if idx == 0 else shoot_dirs[idx - 1]
        print(f"{d:<15} {upstream:<15} {f'{high_val} -> {low_val}':<15}")
    print("-----------------------------------------------")
    print("Done.")


if __name__ == "__main__":
    main()
