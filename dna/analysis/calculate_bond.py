#!/usr/bin/env python3
import argparse
from os import path
from collections import namedtuple
import numpy as np
import pandas as pd
import oxpy
from oxDNA_analysis_tools.UTILS.data_structures import TopInfo, TrajInfo
from oxDNA_analysis_tools.UTILS.oat_multiprocesser import oat_multiprocesser
from oxDNA_analysis_tools.UTILS.RyeReader import describe, get_input_parameter
import time

start_time = time.time()

ComputeContext = namedtuple("ComputeContext", ["traj_info", "top_info", "input_file"])

def compute(ctx: ComputeContext, chunk_size: int, chunk_id: int):
    """Compute total bonds and occupancy per frame for a chunk."""
    with oxpy.Context():
        inp = oxpy.InputFile()
        inp.init_from_filename(ctx.input_file)
        inp["list_type"] = "cells"
        inp["trajectory_file"] = ctx.traj_info.path
        inp["analysis_bytes_to_skip"] = str(ctx.traj_info.idxs[chunk_id*chunk_size].offset)
        inp["confs_to_analyse"] = str(chunk_size)
        inp["analysis_data_output_1"] = """
        {
            name = stdout
            print_every = 1
            col_1 = {
                id = my_obs
                type = hb_list
            }
        }
        """

        backend = oxpy.analysis.AnalysisBackend(inp)

        frame_total_bonds = []
        frame_occupancy = []

        while backend.read_next_configuration():
            pairs = backend.config_info().get_observable_by_id("my_obs").get_output_string(0).strip().split('\n')[1:]
            total_bonds = len(pairs)
            occupancy = total_bonds / (ctx.top_info.nbases / 2) * 100
            frame_total_bonds.append(total_bonds)
            frame_occupancy.append(occupancy)

    return frame_total_bonds, frame_occupancy

def bond_analysis(traj_info: TrajInfo, top_info: TopInfo, inputfile: str, ncpus: int = 1):
    """Run bond analysis across the trajectory, return per-frame totals."""
    ctx = ComputeContext(traj_info, top_info, inputfile)
    all_total_bonds = []
    all_occupancy = []

    def callback(i, r):
        nonlocal all_total_bonds, all_occupancy
        tb, occ = r
        all_total_bonds.extend(tb)
        all_occupancy.extend(occ)

    oat_multiprocesser(traj_info.nconfs, ncpus, compute, callback, ctx)
    return all_total_bonds, all_occupancy

def read_timesteps(traj_file):
    """Read timesteps t = ... from trajectory file."""
    times = []
    with open(traj_file, 'r') as f:
        for line in f:
            if line.startswith("t ="):
                t = float(line.split('=')[1].strip())
                times.append(t)
    return times

def main():
    parser = argparse.ArgumentParser(prog=path.basename(__file__),
                                     description="Calculate total bonds and occupancy % per frame, including wobble pairs")
    parser.add_argument('inputfile', type=str, help="oxDNA input file")
    parser.add_argument('trajectory', type=str, help="Trajectory file")
    parser.add_argument('-p', type=int, default=1, help="Number of cores")
    parser.add_argument('--output', type=str, default="bond_summary.csv", help="Output CSV file name")
    args = parser.parse_args()

    inputfile = args.inputfile
    traj_file = args.trajectory
    ncpus = args.p
    output_file = args.output

    # Read topology and trajectory info
    top_file = get_input_parameter(inputfile, "topology")
    top_info, traj_info = describe(top_file, traj_file)

    # Read actual timesteps from trajectory
    times = read_timesteps(traj_file)

    # Run bond analysis
    total_bonds, occupancy = bond_analysis(traj_info, top_info, inputfile, ncpus)

    # Safety check: ensure number of frames match
    if len(times) != len(total_bonds):
        print(f"WARNING: Number of timesteps ({len(times)}) does not match number of frames ({len(total_bonds)}). Using frame indices as t.")
        times = list(range(len(total_bonds)))

    # Save results to CSV
    df = pd.DataFrame({
        "t": times,
        "total_bonds": total_bonds,
        "occupancy%": occupancy
    })
    df.to_csv(output_file, index=False)
    print(f"Saved bond data to {output_file}")
    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()
