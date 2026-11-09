#!/usr/bin/env python3
import argparse
from os import path
import os
from collections import namedtuple
import re
import operator
import numpy as np
import pandas as pd
import oxpy
from oxDNA_analysis_tools.UTILS.data_structures import TopInfo, TrajInfo
from oxDNA_analysis_tools.UTILS.oat_multiprocesser import (
    oat_multiprocesser,
    get_chunk_size,
)
from oxDNA_analysis_tools.UTILS.RyeReader import describe, get_input_parameter
import time

start_time = time.time()

ComputeContext = namedtuple(
    "ComputeContext",
    ["traj_info", "top_info", "input_file", "indices", "op", "threshold"]
)


def parse_condition(cond_str):
    """
    Parse a condition string like '> 5', '>=10', '== 0' into (op_str, threshold).

    Allowed operators: >, <, >=, <=, ==, !=, =
    '=' is treated as '=='.
    """
    cond_str = cond_str.strip()
    m = re.match(r"(<=|>=|==|!=|<|>|=)\s*(\d+)", cond_str)
    if not m:
        raise ValueError(
            f"Invalid condition '{cond_str}'. "
            "Use forms like '> 5', '>=10', '== 0'."
        )
    op_str, val_str = m.groups()
    if op_str == "=":
        op_str = "=="
    threshold = int(val_str)
    return op_str, threshold


def compute(ctx: ComputeContext, chunk_size: int, chunk_id: int):
    """
    Compute total bonds and occupancy per frame for a chunk,
    restricted to a specific set of nucleotide indices (ctx.indices).

    Additionally: record frames where the total number of bonds satisfies
    the condition (ctx.op, ctx.threshold). The actual extraction of
    configurations is done later in the main process from the raw text
    trajectory.
    """
    indices_set = set(ctx.indices)

    # Map operator string to actual Python comparison function
    op_map = {
        ">": operator.gt,
        "<": operator.lt,
        ">=": operator.ge,
        "<=": operator.le,
        "==": operator.eq,
        "!=": operator.ne,
    }
    if ctx.op not in op_map:
        raise ValueError(f"Unsupported operator '{ctx.op}'")
    compare = op_map[ctx.op]
    threshold = ctx.threshold

    with oxpy.Context():
        inp = oxpy.InputFile()
        inp.init_from_filename(ctx.input_file)

        # Silence oxpy logging
        inp["log_file"] = "/dev/null"

        # Analyze target trajectory and select a chunk
        inp["list_type"] = "cells"
        inp["trajectory_file"] = ctx.traj_info.path
        inp["analysis_bytes_to_skip"] = str(
            ctx.traj_info.idxs[chunk_id * chunk_size].offset
        )
        inp["confs_to_analyse"] = str(chunk_size)

        # Observable setup: hb_list
        inp["analysis_data_output_1"] = (
            "{ \n"
            " name = stdout \n"
            " print_every = 1e10 \n"
            " col_1 = { \n"
            "  id = my_obs \n"
            "  type = hb_list \n"
            " } \n"
            "} "
        )

        backend = oxpy.analysis.AnalysisBackend(inp)

        frame_total_bonds = []
        frame_occupancy = []

        ref_nbases = len(ctx.indices)
        selected_frames = []  # list of (global_frame_idx, step)
        frame_in_chunk = 0

        while backend.read_next_configuration():
            step = backend.config_info().current_step
            raw = (
                backend.config_info()
                .get_observable_by_id("my_obs")
                .get_output_string(step)
            )

            lines = raw.strip().split("\n")
            if not lines:
                total_bonds = 0
                occupancy = 0.0
            else:
                # Skip header line
                pairs_lines = lines[1:]
                total_bonds = 0

                # Only count bonds where both nucleotides are in indices_set
                for p in pairs_lines:
                    p = p.strip()
                    if not p:
                        continue
                    cols = p.split()
                    try:
                        a = int(cols[0])
                        b = int(cols[1])
                    except (ValueError, IndexError):
                        continue
                    if a in indices_set and b in indices_set:
                        total_bonds += 1

                if ref_nbases > 0:
                    occupancy = total_bonds / (ref_nbases / 2.0) * 100.0
                else:
                    occupancy = 0.0

            frame_total_bonds.append(total_bonds)
            frame_occupancy.append(occupancy)

            # Compute global frame index
            global_idx = chunk_id * chunk_size + frame_in_chunk

            # Record frames that satisfy the condition
            if compare(total_bonds, threshold):
                selected_frames.append((global_idx, step))

            frame_in_chunk += 1

    return np.array(frame_total_bonds), np.array(frame_occupancy), selected_frames


def bond_analysis(
    traj_info: TrajInfo,
    top_info: TopInfo,
    inputfile: str,
    indices,
    condition_op: str,
    condition_threshold: int,
    ncpus: int = 1,
):
    """
    Run bond analysis across the trajectory for a specific set of indices.

    Additionally, collect the list of (frame_idx, step) for frames
    satisfying the condition (condition_op, condition_threshold).
    """
    ctx = ComputeContext(
        traj_info, top_info, inputfile, indices, condition_op, condition_threshold
    )

    total_bonds = np.empty(traj_info.nconfs)
    occupancy = np.empty(traj_info.nconfs)

    chunk_size = get_chunk_size()
    selected_frames_all = []  # (frame_idx, step)

    def callback(i, r):
        nonlocal total_bonds, occupancy, selected_frames_all
        tb, occ, selected_chunk = r
        start = i * chunk_size
        end = start + len(tb)
        total_bonds[start:end] = tb
        occupancy[start:end] = occ
        selected_frames_all.extend(selected_chunk)

    oat_multiprocesser(traj_info.nconfs, ncpus, compute, callback, ctx)
    return total_bonds, occupancy, selected_frames_all


def read_timesteps(traj_file):
    """Read timesteps t = ... from trajectory file."""
    times = []
    with open(traj_file, "r") as f:
        for line in f:
            if line.startswith("t ="):
                t = float(line.split("=", 1)[1].strip().split()[0])
                times.append(t)
    return times


def read_indices_file(idx_file):
    """Read comma-separated nucleotide indices from a file (0-based)."""
    with open(idx_file, "r") as f:
        content = f.read()
    tokens = content.replace("\n", ",").split(",")
    indices = [int(tok.strip()) for tok in tokens if tok.strip() != ""]
    if len(indices) == 0:
        raise ValueError(f"No valid indices found in {idx_file}")
    return indices


def write_configs_from_text(traj_file, selected_steps, output_dir):
    """
    Parse the plain-text trajectory file and, for every block starting with
    't = <step>', if <step> is in selected_steps, write that block to
    '<output_dir>/conf_<step>.dat'.

    Returns the number of distinct configurations actually written.
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    if not selected_steps:
        print("No steps matched the condition; no configurations written.")
        return 0

    # Normalize all selected steps to canonical integer strings
    selected_strs = {str(int(s)) for s in selected_steps if s is not None}

    current_lines = []
    current_step_str = None
    written_steps = set()

    def flush_current():
        if current_step_str is None:
            return
        if current_step_str in selected_strs and current_lines:
            out_name = f"conf_{current_step_str}.dat"
            out_path = os.path.join(output_dir, out_name)
            with open(out_path, "w") as out:
                out.writelines(current_lines)
            written_steps.add(current_step_str)

    with open(traj_file, "r") as f:
        for line in f:
            if line.startswith("t ="):
                # End of previous block
                flush_current()

                # Start new block
                current_lines = [line]
                try:
                    # Take the first token after '='; e.g. "6000000"
                    step_token = line.split("=", 1)[1].strip().split()[0]
                    # Normalize as integer string
                    current_step_str = str(int(float(step_token)))
                except Exception:
                    current_step_str = None
                continue

            # Lines belonging to the current configuration
            current_lines.append(line)

        # Flush last block at EOF
        flush_current()

    print(f"Wrote configurations for {len(written_steps)} distinct steps into '{output_dir}'.")
    return len(written_steps)


def main():
    parser = argparse.ArgumentParser(
        prog=path.basename(__file__),
        description=(
            "Calculate total bonds and occupancy %% per frame for a specific set of nucleotides,\n"
            "including wobble pairs, and save configurations for frames\n"
            "where the HB count matches a given condition. Indices are provided as a CSV file."
        ),
    )
    parser.add_argument("inputfile", type=str, help="oxDNA input file")
    parser.add_argument("trajectory", type=str, help="Trajectory file")
    parser.add_argument(
        "indexes",
        type=str,
        help=(
            "Path to a CSV-like file with nucleotide indices "
            "(comma-separated, e.g. '1,3,6,4,2'). "
            "Assumed 0-based; only bonds between these indices are counted."
        ),
    )
    parser.add_argument(
        "condition",
        type=str,
        help=(
            "Condition on the HB count for saving configurations, e.g. "
            "'> 5', '>=10', '== 0'. If satisfied, conf_<step>.dat is written."
        ),
    )
    parser.add_argument("-p", type=int, default=1, help="Number of cores")
    parser.add_argument(
        "--output",
        type=str,
        default="bond_summary.csv",
        help="Output CSV file name",
    )
    args = parser.parse_args()

    inputfile = args.inputfile
    traj_file = args.trajectory
    idx_file = args.indexes
    cond_str = args.condition
    ncpus = args.p
    output_file = args.output

    # Parse condition
    cond_op, cond_thr = parse_condition(cond_str)

    # Build output directory name, e.g. HB_KL>5 from condition '> 5'
    # Remove spaces to keep a clean directory name.
    cond_for_name = cond_str.replace(" ", "")
    output_dir = f"HB_KL{cond_for_name}"
    os.makedirs(output_dir, exist_ok=True)

    # Read indices (required)
    indices = read_indices_file(idx_file)

    # Read topology and trajectory info
    top_file = get_input_parameter(inputfile, "topology")
    top_info, traj_info = describe(top_file, traj_file)

    # Read actual timesteps from trajectory
    times = read_timesteps(traj_file)

    # Run bond analysis (also collecting frames that match the condition)
    total_bonds, occupancy, selected_frames = bond_analysis(
        traj_info,
        top_info,
        inputfile,
        indices=indices,
        condition_op=cond_op,
        condition_threshold=cond_thr,
        ncpus=ncpus,
    )

    # Extract the step values from (frame_idx, step)
    selected_steps = [step for _, step in selected_frames]
    distinct_selected_steps = set(selected_steps)

    total_frames = len(total_bonds)  # number of HB calculations
    satisfied_count = len(distinct_selected_steps)
    probability = (satisfied_count / total_frames) if total_frames > 0 else 0.0

    print(f"Condition: {cond_str}")
    print(f"Total frames (HB calculations): {total_frames}")
    print(f"Distinct timesteps satisfying condition: {satisfied_count}")
    print(f"Probability: {probability:.6f} ({satisfied_count}/{total_frames})")

    # Write plain-text configs for those steps into output_dir
    _ = write_configs_from_text(traj_file, selected_steps, output_dir)

    # Safety check: ensure number of frames match
    if len(times) != len(total_bonds):
        print(
            f"WARNING: Number of timesteps ({len(times)}) does not match number "
            f"of frames ({len(total_bonds)}). Using frame indices as t."
        )
        times = list(range(len(total_bonds)))

    # Save results to CSV (you can remove this block if you only care about conf_* files)
    df = pd.DataFrame(
        {
            "t": times,
            "total_bonds": total_bonds.astype(int),
            "occupancy%": occupancy,
        }
    )
    df.to_csv(output_file, index=False)
    print(f"Saved bond data to {output_file}")

    # Save summary log in the same directory as conf files
    log_path = os.path.join(output_dir, "hb_conf.log")
    with open(log_path, "w") as logf:
        logf.write(f"Condition: {cond_str}\n")
        logf.write(f"Total_frames(HB_calculations): {total_frames}\n")
        logf.write(f"Distinct_timesteps_satisfying_condition: {satisfied_count}\n")
        logf.write(f"Probability: {probability:.6f} ({satisfied_count}/{total_frames})\n")
        logf.write(f"Config_directory: {output_dir}\n")

    print(f"Wrote summary log to {log_path}")
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
