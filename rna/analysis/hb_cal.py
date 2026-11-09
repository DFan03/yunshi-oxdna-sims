#!/usr/bin/env python3
import argparse
from os import path
from collections import namedtuple
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
    "ComputeContext", ["traj_info", "top_info", "input_file", "indices"]
)


def compute(ctx: ComputeContext, chunk_size: int, chunk_id: int):
    """Compute total bonds and occupancy per frame for a chunk,
    restricted to a specific set of nucleotide indices (ctx.indices)."""
    indices_set = set(ctx.indices)

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

        # Same observable setup as in the package example
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

        while backend.read_next_configuration():
            step = backend.config_info().current_step
            raw = (
                backend.config_info()
                .get_observable_by_id("my_obs")
                .get_output_string(step)
            )

            lines = raw.strip().split("\n")
            if not lines:
                frame_total_bonds.append(0)
                frame_occupancy.append(0.0)
                continue

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

    return np.array(frame_total_bonds), np.array(frame_occupancy)


def bond_analysis(
    traj_info: TrajInfo,
    top_info: TopInfo,
    inputfile: str,
    indices,
    ncpus: int = 1,
):
    """Run bond analysis across the trajectory for a specific set of indices."""
    ctx = ComputeContext(traj_info, top_info, inputfile, indices)

    total_bonds = np.empty(traj_info.nconfs)
    occupancy = np.empty(traj_info.nconfs)

    chunk_size = get_chunk_size()

    def callback(i, r):
        nonlocal total_bonds, occupancy
        tb, occ = r
        start = i * chunk_size
        end = start + len(tb)
        total_bonds[start:end] = tb
        occupancy[start:end] = occ

    oat_multiprocesser(traj_info.nconfs, ncpus, compute, callback, ctx)
    return total_bonds, occupancy


def read_timesteps(traj_file):
    """Read timesteps t = ... from trajectory file."""
    times = []
    with open(traj_file, "r") as f:
        for line in f:
            if line.startswith("t ="):
                t = float(line.split("=")[1].strip())
                times.append(t)
    return times


def read_indices_file(idx_file):
    """Read comma-separated nucleotide indices from a file.

    The file may contain commas and/or newlines, e.g.:
      1,3,6,4,2

    Assumes indices are 0-based, matching hb_list.
    """
    with open(idx_file, "r") as f:
        content = f.read()
    tokens = content.replace("\n", ",").split(",")
    indices = [int(tok.strip()) for tok in tokens if tok.strip() != ""]
    if len(indices) == 0:
        raise ValueError(f"No valid indices found in {idx_file}")
    return indices


def main():
    parser = argparse.ArgumentParser(
        prog=path.basename(__file__),
        description=(
            "Calculate total bonds and occupancy %% per frame for a specific set of nucleotides,\n"
            "including wobble pairs. Indices are provided as a CSV file."
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
    ncpus = args.p
    output_file = args.output

    # Read indices (required)
    indices = read_indices_file(idx_file)

    # Read topology and trajectory info
    top_file = get_input_parameter(inputfile, "topology")
    top_info, traj_info = describe(top_file, traj_file)

    # Read actual timesteps from trajectory
    times = read_timesteps(traj_file)

    # Run bond analysis
    total_bonds, occupancy = bond_analysis(
        traj_info, top_info, inputfile, indices=indices, ncpus=ncpus
    )

    # Safety check: ensure number of frames match
    if len(times) != len(total_bonds):
        print(
            f"WARNING: Number of timesteps ({len(times)}) does not match number "
            f"of frames ({len(total_bonds)}). Using frame indices as t."
        )
        times = list(range(len(total_bonds)))

    # Save results to CSV
    df = pd.DataFrame(
        {
            "t": times,
            "total_bonds": total_bonds.astype(int),
            "occupancy%": occupancy,
        }
    )
    df.to_csv(output_file, index=False)
    print(f"Saved bond data to {output_file}")
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
