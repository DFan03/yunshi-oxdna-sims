"""
setup_umbrella.py
-----------------
Reads umbrella_params.txt and sets up all umbrella sampling windows.

Directory structure created:
  umbrella/
    window_00/
      equil/
        input.txt
        forces.txt
      prod/
        input.txt
        forces.txt
    window_01/
      ...
  run_all.sh   ← sequential run script

Usage:
    python setup_umbrella.py umbrella_params.txt
"""

import os
import sys
import csv
import shutil
import numpy as np


# ── Helpers ──────────────────────────────────────────────────────────────────

def parse_params(filepath):
    """Parse key = value param file, ignoring # comments."""
    params = {}
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if '=' not in line:
                continue
            key, _, val = line.partition('=')
            params[key.strip()] = val.strip()
    return params


def read_csv_indices(filepath):
    """Read integer nucleotide indices from a CSV (0-based, one per line/cell)."""
    indices = []
    with open(filepath, newline='') as f:
        reader = csv.reader(f)
        for row in reader:
            for val in row:
                val = val.strip()
                if val and not val.startswith('#'):
                    indices.append(int(val))
    return indices


def read_template(filepath):
    """Read the oxDNA template input file as a string."""
    with open(filepath) as f:
        return f.read()


def write_input(template, replacements, out_path):
    """
    Write an oxDNA input file by replacing keys in the template.
    replacements: dict of {key: new_value_string}
    Keys already present in the template are overwritten;
    keys not present are appended.
    """
    lines = template.splitlines()
    used_keys = set()
    new_lines = []

    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith('#'):
            new_lines.append(line)
            continue
        if '=' in line:
            key = line.split('=')[0].strip()
            if key in replacements:
                new_lines.append(f"{key} = {replacements[key]}")
                used_keys.add(key)
                continue
        new_lines.append(line)

    # Append any keys that weren't in the template
    for key, val in replacements.items():
        if key not in used_keys:
            new_lines.append(f"{key} = {val}")

    with open(out_path, 'w') as f:
        f.write('\n'.join(new_lines) + '\n')


def write_forces(group_a, group_b, stiff, r0, out_path):
    """Write oxDNA com_force external forces file for one window."""
    a_str = ' '.join(str(i) for i in group_a)
    b_str = ' '.join(str(i) for i in group_b)
    content = (
        "{\n"
        "type = com\n"
        f"group_particle_1 = {a_str}\n"
        f"group_particle_2 = {b_str}\n"
        f"stiff = {stiff}\n"
        f"r0 = {r0:.6f}\n"
        "PBC = 1\n"
        "}\n"
    )
    with open(out_path, 'w') as f:
        f.write(content)


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    if len(sys.argv) != 2:
        print("Usage: python setup_umbrella.py umbrella_params.txt")
        sys.exit(1)

    param_file = sys.argv[1]
    base_dir = os.path.dirname(os.path.abspath(param_file))
    p = parse_params(param_file)

    # ── Load parameters ──
    topology    = os.path.join(base_dir, p['topology'])
    conf_file   = os.path.join(base_dir, p['conf_file'])
    group_a_csv = os.path.join(base_dir, p['group_a_csv'])
    group_b_csv = os.path.join(base_dir, p['group_b_csv'])
    template_f  = os.path.join(base_dir, p['md_template'])
    oxdna_bin   = p['oxdna_bin']

    xi_min      = float(p['xi_min'])
    xi_max      = float(p['xi_max'])
    n_windows   = int(p['n_windows'])
    stiff       = float(p['stiff'])

    equil_steps = int(p['equil_steps'])
    prod_steps  = int(p['prod_steps'])

    temperature        = p.get('temperature', '20C')
    salt_concentration = p.get('salt_concentration', '1')
    seq_dep_file       = p.get('seq_dep_file', '')

    # ── Load nucleotide indices ──
    group_a = read_csv_indices(group_a_csv)
    group_b = read_csv_indices(group_b_csv)
    print(f"Group A: {len(group_a)} nucleotides")
    print(f"Group B: {len(group_b)} nucleotides")

    # ── Window target distances ──
    xi_values = np.linspace(xi_min, xi_max, n_windows)
    print(f"Windows: {n_windows}  |  xi range: {xi_min} – {xi_max} (spacing {xi_values[1]-xi_values[0]:.4f})")

    # ── Load template ──
    template = read_template(template_f)

    # ── Output root ──
    out_root = os.path.join(base_dir, 'umbrella')
    os.makedirs(out_root, exist_ok=True)

    # Copy topology to umbrella root for easy reference
    shutil.copy(topology, os.path.join(out_root, os.path.basename(topology)))

    run_commands = []

    for i, xi0 in enumerate(xi_values):
        win_name = f"window_{i:02d}"
        win_dir  = os.path.join(out_root, win_name)

        equil_dir = os.path.join(win_dir, 'equil')
        prod_dir  = os.path.join(win_dir, 'prod')
        os.makedirs(equil_dir, exist_ok=True)
        os.makedirs(prod_dir,  exist_ok=True)

        # ── Forces file (same xi0 for both stages) ──
        write_forces(group_a, group_b, stiff, xi0,
                     os.path.join(equil_dir, 'forces.txt'))
        write_forces(group_a, group_b, stiff, xi0,
                     os.path.join(prod_dir,  'forces.txt'))

        # Relative topology path (both stages are one level below win_dir)
        rel_topology = os.path.relpath(topology, equil_dir)

        # ── Equilibration input ──
        equil_conf = os.path.relpath(conf_file, equil_dir)

        # print intervals: save ~1000 frames
        equil_print_conf  = max(1, equil_steps // 1000)
        equil_print_energy = max(1, equil_steps // 10000)

        equil_replacements = {
            'topology':               rel_topology,
            'conf_file':              equil_conf,
            'T':                      temperature,
            'salt_concentration':     salt_concentration,
            'seq_dep_file':           seq_dep_file,
            'steps':                  str(equil_steps),
            'print_conf_interval':    str(equil_print_conf),
            'print_energy_every':     str(equil_print_energy),
            'external_forces':        '1',
            'external_forces_file':   'forces.txt',
            'lastconf_file':          'last_conf.dat',
            'trajectory_file':        'trajectory.dat',
            'energy_file':            'energy.dat',
            'refresh_vel':            '1',
            'restart_step_counter':   '1',
        }
        write_input(template, equil_replacements,
                    os.path.join(equil_dir, 'input.txt'))

        # ── Production input ──
        # Starts from equilibration last_conf
        prod_conf = os.path.relpath(
            os.path.join(equil_dir, 'last_conf.dat'), prod_dir)
        rel_topology_prod = os.path.relpath(topology, prod_dir)

        prod_print_conf   = max(1, prod_steps // 2000)
        prod_print_energy = max(1, prod_steps // 10000)

        prod_replacements = {
            'topology':               rel_topology_prod,
            'conf_file':              prod_conf,
            'T':                      temperature,
            'salt_concentration':     salt_concentration,
            'seq_dep_file':           seq_dep_file,
            'steps':                  str(prod_steps),
            'print_conf_interval':    str(prod_print_conf),
            'print_energy_every':     str(prod_print_energy),
            'external_forces':        '1',
            'external_forces_file':   'forces.txt',
            'lastconf_file':          'last_conf.dat',
            'trajectory_file':        'trajectory.dat',
            'energy_file':            'energy.dat',
            'refresh_vel':            '0',
            'restart_step_counter':   '1',
        }
        write_input(template, prod_replacements,
                    os.path.join(prod_dir, 'input.txt'))

        # ── Run commands ──
        run_commands.append(f"# Window {i:02d}  xi0 = {xi0:.4f}")
        run_commands.append(f"echo '[{win_name}] Starting equilibration...'")
        run_commands.append(f"cd {equil_dir} && {oxdna_bin} input.txt > run.log 2>&1")
        run_commands.append(f"echo '[{win_name}] Starting production...'")
        run_commands.append(f"cd {prod_dir}  && {oxdna_bin} input.txt > run.log 2>&1")
        run_commands.append("")

        print(f"  [{win_name}]  xi0 = {xi0:.4f}")

    # ── Write run_all.sh ──
    run_script = os.path.join(out_root, 'run_all.sh')
    with open(run_script, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write("# Auto-generated by setup_umbrella.py\n")
        f.write("# Runs all umbrella windows sequentially (equil then prod per window)\n\n")
        f.write("set -e  # stop on first error\n\n")
        f.write('\n'.join(run_commands) + '\n')
        f.write('\necho "All windows complete."\n')

    os.chmod(run_script, 0o755)

    print(f"\nSetup complete.")
    print(f"  Output directory : {out_root}")
    print(f"  Run script       : {run_script}")
    print(f"  Launch with      : bash {run_script}")


if __name__ == '__main__':
    main()

