"""
persistence_length_csv.py - SINGLE HELIX MODE

Treats the entire structure as ONE continuous helix and reports ONE
persistence length value. All segments are processed together with a
single global axis.
"""

from os import path
from typing import Dict, List, Tuple
from collections import namedtuple
import argparse
import numpy as np
import matplotlib.pyplot as plt
import oxpy
from oxDNA_analysis_tools.UTILS.RyeReader import describe
from oxDNA_analysis_tools.UTILS.oat_multiprocesser import oat_multiprocesser
from oxDNA_analysis_tools.UTILS.data_structures import TrajInfo
from oxDNA_analysis_tools.UTILS.logger import log, logger_settings
from oxDNA_analysis_tools.UTILS.constants import FIG_DPI

import time
start_time = time.time()

ComputeContext = namedtuple("ComputeContext", [
    "traj_info",
    "input_file",
    "nuc_indices",
    "pair_dict",
    "valid_pairs",
    "step_signs",
])

BP_DISTANCE_CUTOFF = 0.65


def load_indices_from_csv(csv_path: str) -> List[int]:
    indices = []
    with open(csv_path, newline='') as f:
        for line in f:
            for token in line.replace(',', ' ').split():
                try:
                    indices.append(int(token))
                except ValueError:
                    continue
    if len(indices) < 2:
        raise ValueError(f"CSV file '{csv_path}' must contain at least 2 nucleotide indices.")
    return indices


def detect_segments(nuc_indices: List[int], gap_threshold: int = 2) -> List[Tuple[int, int]]:
    segments, seg_start = [], 0
    for i in range(1, len(nuc_indices)):
        if abs(nuc_indices[i] - nuc_indices[i - 1]) > gap_threshold:
            segments.append((seg_start, i - 1))
            seg_start = i
    segments.append((seg_start, len(nuc_indices) - 1))
    return segments


def make_pos_to_seg(nuc_indices, segments):
    pos_to_seg = {}
    for seg_idx, (s, e) in enumerate(segments):
        for p in range(s, e + 1):
            pos_to_seg[p] = seg_idx
    return pos_to_seg


def read_hb_map(inp_file: str, traj_info: TrajInfo, frame_idx: int) -> Dict[int, int]:
    offset = str(traj_info.idxs[frame_idx].offset)

    def _pass(overrides):
        result = {}
        with oxpy.Context():
            inp = oxpy.InputFile()
            inp.init_from_filename(inp_file)
            inp["log_file"] = "/dev/null"
            inp["list_type"] = "cells"
            inp["trajectory_file"] = traj_info.path
            inp["analysis_bytes_to_skip"] = offset
            inp["confs_to_analyse"] = "1"
            inp["analysis_data_output_1"] = (
                '{ \n name = stdout \n print_every = 1e10 \n'
                ' col_1 = { \n id = my_obs \n type = hb_list \n } \n }'
            )
            for k, v in overrides.items():
                inp[k] = v
            backend = oxpy.analysis.AnalysisBackend(inp)
            backend.read_next_configuration()
            raw = (
                backend.config_info()
                .get_observable_by_id("my_obs")
                .get_output_string(backend.config_info().current_step)
                .strip().split('\n')[1:]
            )
        for p in raw:
            fields = p.split()
            if len(fields) < 2:
                continue
            p1, p2 = int(fields[0]), int(fields[1])
            result[p1] = p2
            result[p2] = p1
        return result

    hb_map = {}
    hb_map.update(_pass({}))
    hb_map.update(_pass({"use_average_seq": "0", "interaction_type": "RNA"}))
    return hb_map


def expected_partner_positions(p, confirmed_anchors, n):
    candidates = []
    for (a, b) in confirmed_anchors:
        for cand in (b - (p - a), a - (p - b)):
            if 0 <= cand < n and cand != p:
                candidates.append(cand)
    return candidates


def find_pairs_by_geometry_register(inp_file, traj_info, frame_idx, nuc_indices, already_paired, confirmed_anchors):
    offset = str(traj_info.idxs[frame_idx].offset)
    n = len(nuc_indices)
    new_pairs = {}

    unresolved = [nuc for nuc in nuc_indices if nuc not in already_paired]
    if not unresolved:
        return new_pairs

    with oxpy.Context():
        inp = oxpy.InputFile()
        inp.init_from_filename(inp_file)
        inp["log_file"] = "/dev/null"
        inp["list_type"] = "cells"
        inp["trajectory_file"] = traj_info.path
        inp["analysis_bytes_to_skip"] = offset
        inp["confs_to_analyse"] = "1"
        inp["analysis_data_output_1"] = (
            '{ \n name = stdout \n print_every = 1e10 \n'
            ' col_1 = { \n id = my_step \n type = step \n } \n }'
        )
        backend = oxpy.analysis.AnalysisBackend(inp)
        backend.read_next_configuration()
        conf = backend.config_info()
        box = np.array(conf.box.box_sides)
        positions = {nuc: np.array(conf.particles()[nuc].base_site()) for nuc in nuc_indices}

    idx_position = {nuc: pos for pos, nuc in enumerate(nuc_indices)}
    committed = {}

    for u in unresolved:
        if u in committed:
            continue
        p = idx_position[u]
        pos_u = positions[u]

        cand_positions = (expected_partner_positions(p, confirmed_anchors, n)
                          if confirmed_anchors else
                          [idx_position[v] for v in nuc_indices
                           if v != u and v not in already_paired and v not in committed])

        best_dist, best_partner = BP_DISTANCE_CUTOFF, None
        for cand_pos in cand_positions:
            v = nuc_indices[cand_pos]
            if v in already_paired or v in committed:
                continue
            diff = pos_u - positions[v]
            diff -= box * np.rint(diff / box)
            dist = np.linalg.norm(diff)
            if dist < best_dist:
                best_dist, best_partner = dist, v

        if best_partner is not None:
            new_pairs[u] = best_partner
            new_pairs[best_partner] = u
            committed[u] = best_partner
            committed[best_partner] = u

    return new_pairs


def build_pair_dict(inp_file, traj_info, nuc_indices):
    idx_position = {nuc: pos for pos, nuc in enumerate(nuc_indices)}
    n = len(nuc_indices)

    mid_frame = traj_info.nconfs // 2
    log(f"Building pair_dict from frame {mid_frame} (middle of {traj_info.nconfs} frames).")

    hb_map = read_hb_map(inp_file, traj_info, mid_frame)

    pair_dict = {}
    for nuc in nuc_indices:
        if nuc in hb_map:
            partner = hb_map[nuc]
            if partner in idx_position:
                pair_dict[nuc] = partner
                pair_dict[partner] = nuc
    log(f"Stage 1 (hb_list): {len(pair_dict) // 2:3d} pairs confirmed.")

    confirmed_anchors = list({
        (min(idx_position[nuc], idx_position[pair_dict[nuc]]),
         max(idx_position[nuc], idx_position[pair_dict[nuc]]))
        for nuc in pair_dict
        if nuc in idx_position and pair_dict[nuc] in idx_position
    })

    geo_pairs = find_pairs_by_geometry_register(inp_file, traj_info, mid_frame, nuc_indices,
                                                set(pair_dict.keys()), confirmed_anchors)
    pair_dict.update(geo_pairs)
    log(f"Stage 2 (geometry): {len(pair_dict) // 2:3d} pairs confirmed (+{len(geo_pairs) // 2} from geometry).")

    confirmed_anchors = list({
        (min(idx_position[nuc], idx_position[pair_dict[nuc]]),
         max(idx_position[nuc], idx_position[pair_dict[nuc]]))
        for nuc in pair_dict
        if nuc in idx_position and pair_dict[nuc] in idx_position
    })

    segments = detect_segments(nuc_indices)
    pos_to_seg = make_pos_to_seg(nuc_indices, segments)
    log(f"Segments: {[(nuc_indices[s], nuc_indices[e]) for s, e in segments]}")

    unresolved = [i for i in range(n) if nuc_indices[i] not in pair_dict]
    for p in unresolved:
        nuc_p = nuc_indices[p]
        seg_p = pos_to_seg[p]
        sorted_anchors = sorted(confirmed_anchors,
            key=lambda ab: (0 if pos_to_seg[ab[0]] == seg_p or pos_to_seg[ab[1]] == seg_p else 1))
        resolved = False
        for (a, b) in sorted_anchors:
            for cand_pos in (b - (p - a), a - (p - b)):
                if cand_pos < 0 or cand_pos >= n:
                    continue
                cand_nuc = nuc_indices[cand_pos]
                if cand_nuc in pair_dict:
                    continue
                if hb_map.get(nuc_p) == cand_nuc:
                    pair_dict[nuc_p] = cand_nuc
                    pair_dict[cand_nuc] = nuc_p
                    resolved = True
                    break
            if resolved:
                break
    log(f"Stage 3 (anchor infer): {len(pair_dict) // 2:3d} pairs confirmed.")

    paired = sum(1 for nuc in nuc_indices if nuc in pair_dict)
    unpaired = n - paired
    log(f"Final: {paired} paired, {unpaired} truly unpaired.")
    if unpaired > 0:
        log(f"Truly unpaired IDs: {[nuc for nuc in nuc_indices if nuc not in pair_dict]}")

    return pair_dict


def get_valid_pairs_and_signs(inp_file: str, traj_info: TrajInfo, nuc_indices: List[int], pair_dict: Dict[int, int]) -> Tuple[List[Tuple[int, int]], np.ndarray]:
    """Build ALL step pairs and compute a SINGLE global helix axis."""
    paired_in_order = [nuc for nuc in nuc_indices if nuc in pair_dict]

    if len(paired_in_order) < 2:
        raise RuntimeError("Fewer than 2 paired nucleotides. Cannot compute PL.")

    valid_pairs = [(paired_in_order[i], paired_in_order[i + 1]) for i in range(len(paired_in_order) - 1)]
    log(f"Valid step pairs (all, including cross-segment): {len(valid_pairs)}")

    mid_frame = traj_info.nconfs // 2
    offset = str(traj_info.idxs[mid_frame].offset)

    with oxpy.Context():
        inp = oxpy.InputFile()
        inp.init_from_filename(inp_file)
        inp["log_file"] = "/dev/null"
        inp["list_type"] = "cells"
        inp["trajectory_file"] = traj_info.path
        inp["analysis_bytes_to_skip"] = offset
        inp["confs_to_analyse"] = "1"
        inp["analysis_data_output_1"] = (
            '{ \n name = stdout \n print_every = 1e10 \n'
            ' col_1 = { \n id = my_step \n type = step \n } \n }'
        )
        backend = oxpy.analysis.AnalysisBackend(inp)
        backend.read_next_configuration()
        conf = backend.config_info()
        box = np.array(conf.box.box_sides)

        raw_vecs = np.zeros((len(valid_pairs), 3))
        for v, (nuc_a, nuc_b) in enumerate(valid_pairs):
            fA = np.array(conf.particles()[nuc_a].base_site())
            fB = np.array(conf.particles()[pair_dict[nuc_a]].base_site())
            sA = np.array(conf.particles()[nuc_b].base_site())
            sB = np.array(conf.particles()[pair_dict[nuc_b]].base_site())
            r = (sA + sB) / 2 - (fA + fB) / 2
            r -= box * np.rint(r / box)
            norm = np.linalg.norm(r)
            if norm > 1e-10:
                raw_vecs[v] = r / norm

    axis = np.mean(raw_vecs, axis=0)
    axis_norm = np.linalg.norm(axis)
    if axis_norm > 1e-10:
        axis /= axis_norm

    log(f"Global helix axis (mean): [{axis[0]:.3f}, {axis[1]:.3f}, {axis[2]:.3f}]")

    step_signs = np.where(raw_vecs @ axis >= 0, 1.0, -1.0)
    n_flipped = int((step_signs < 0).sum())
    log(f"Steps flipped for axis alignment: {n_flipped} / {len(valid_pairs)}")

    return valid_pairs, step_signs


def get_r_step(conf, nuc_a: int, nuc_b: int, pair_dict: Dict) -> np.ndarray:
    box = np.array(conf.box.box_sides)
    fA = np.array(conf.particles()[nuc_a].base_site())
    fB = np.array(conf.particles()[pair_dict[nuc_a]].base_site())
    sA = np.array(conf.particles()[nuc_b].base_site())
    sB = np.array(conf.particles()[pair_dict[nuc_b]].base_site())
    r = (sA + sB) / 2 - (fA + fB) / 2
    r -= box * np.rint(r / box)
    return r


def compute(ctx: ComputeContext, chunk_size: int, chunk_id: int):
    pair_dict = ctx.pair_dict
    valid_pairs = ctx.valid_pairs
    step_signs = ctx.step_signs
    n_valid = len(valid_pairs)

    l0 = 0.0
    correlations = np.zeros(n_valid)
    correlations_counter = np.zeros(n_valid)

    with oxpy.Context():
        inp = oxpy.InputFile()
        inp.init_from_filename(ctx.input_file)
        inp["log_file"] = "/dev/null"
        inp["list_type"] = "cells"
        inp["trajectory_file"] = ctx.traj_info.path
        inp["analysis_bytes_to_skip"] = str(ctx.traj_info.idxs[chunk_id * chunk_size].offset)
        inp["confs_to_analyse"] = str(chunk_size)
        inp["analysis_data_output_1"] = (
            '{ \n name = stdout \n print_every = 1e10 \n'
            ' col_1 = { \n id = my_step \n type = step \n } \n }'
        )

        backend = oxpy.analysis.AnalysisBackend(inp)

        while backend.read_next_configuration():
            conf = backend.config_info()

            step_vecs = np.empty((n_valid, 3))
            step_valid = np.ones(n_valid, dtype=bool)

            for v, (nuc_a, nuc_b) in enumerate(valid_pairs):
                r = get_r_step(conf, nuc_a, nuc_b, pair_dict)
                norm = np.linalg.norm(r)
                if norm < 1e-10:
                    step_valid[v] = False
                    continue
                l0 += norm
                step_vecs[v] = step_signs[v] * (r / norm)

            sv_idx = np.where(step_valid)[0]
            if len(sv_idx) < 2:
                continue
            
            sv = step_vecs[sv_idx]
            nv = sv.shape[0]
            dots = sv @ sv.T

            for offset in range(nv):
                correlations[offset] += np.trace(dots, offset=offset)
                correlations_counter[offset] += nv - offset

    return (l0, correlations, correlations_counter)


def persistence_length(traj_info, inp_file, nuc_indices, ncpus=1):
    pair_dict = build_pair_dict(inp_file, traj_info, nuc_indices)
    valid_pairs, step_signs = get_valid_pairs_and_signs(inp_file, traj_info, nuc_indices, pair_dict)

    n_valid = len(valid_pairs)
    ctx = ComputeContext(traj_info, inp_file, nuc_indices, pair_dict, valid_pairs, step_signs)

    l0 = 0.0
    correlations = np.zeros(n_valid)
    correlations_counter = np.zeros(n_valid)

    def callback(i, r):
        nonlocal l0, correlations, correlations_counter
        l0 += r[0]
        correlations += r[1]
        correlations_counter += r[2]

    oat_multiprocesser(traj_info.nconfs, ncpus, compute, callback, ctx)

    l0 /= (traj_info.nconfs * n_valid)
    
    with np.errstate(invalid='ignore'):
        correlations = np.where(correlations_counter > 0,
                                correlations / correlations_counter,
                                np.nan)
    return (l0, correlations)


def fit_PL(correlations: np.ndarray, plt_name: str) -> float:
    x_all = np.arange(len(correlations))
    usable = np.isfinite(correlations) & (correlations > 0)

    if not usable[0]:
        raise RuntimeError("Correlation at offset 0 is not positive.")

    non_usable = np.where(~usable)[0]
    cutoff = non_usable[0] if len(non_usable) > 0 else len(correlations)

    if cutoff < 2:
        raise RuntimeError(f"Only {cutoff} usable correlation value(s).")

    x = x_all[:cutoff]
    log_corr = np.log(correlations[:cutoff])
    A, B = np.polyfit(x, log_corr, 1)
    pl = -1.0 / A

    log(f"Fit used offsets 0–{cutoff-1} ({cutoff} points out of {len(correlations)})")

    safe_log = np.where(usable, np.log(np.where(usable, correlations, np.nan)), np.nan)
    fig, ax = plt.subplots()
    ax.scatter(x_all, safe_log, alpha=0.6, label='ln(correlation)')
    if (~usable).any():
        ax.scatter(x_all[~usable], np.full((~usable).sum(), np.nanmin(safe_log)),
                   marker='x', color='red', alpha=0.5, label='non-positive/NaN (excluded)')
    trend = np.poly1d([A, B])
    ax.plot(x, trend(x), color='tab:orange', label=f'Fit y = (-x/{pl:.1f})+{B:.3f}')
    ax.spines[['top', 'right']].set_visible(False)
    ax.set_xlabel('Offset (base-pair steps along helix)')
    ax.set_ylabel('ln(correlation)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(plt_name, dpi=FIG_DPI)
    plt.close(fig)
    log(f"Saving figure to {plt_name}")
    return pl


def cli_parser(prog="persistence_length_csv.py"):
    parser = argparse.ArgumentParser(prog=prog, description=(
        "Calculates persistence length treating the entire helix as ONE continuous structure."
    ))
    parser.add_argument('traj_file', type=str, nargs=1)
    parser.add_argument('input', type=str, nargs=1)
    parser.add_argument('csv_file', type=str, nargs=1)
    parser.add_argument('-p', '--parallel', metavar='num_cpus', nargs=1, type=int, dest='parallel')
    parser.add_argument('-d', '--data', metavar='data_file', nargs=1)
    parser.add_argument('-n', '--plot_name', nargs=1)
    parser.add_argument('-q', '--quiet', dest='quiet', action='store_const', const=True, default=False)
    return parser


def main():
    parser = cli_parser(path.basename(__file__))
    args = parser.parse_args()

    logger_settings.set_quiet(args.quiet)
    from oxDNA_analysis_tools.config import check
    check(["python", "matplotlib", "numpy", "oxpy"])

    traj_file = args.traj_file[0]
    inp_file = args.input[0]
    csv_file = args.csv_file[0]

    nuc_indices = load_indices_from_csv(csv_file)
    log(f"Loaded {len(nuc_indices)} nucleotide indices from '{csv_file}'")
    log(f"First: {nuc_indices[0]}, Last: {nuc_indices[-1]}")

    ncpus = args.parallel[0] if args.parallel else 1

    ti, di = describe(None, traj_file)
    l0, correlations = persistence_length(di, inp_file, nuc_indices, ncpus)

    if args.data:
        log(f"Writing correlations to {args.data[0]}")
        with open(args.data[0], 'w+') as f:
            for i, c in enumerate(correlations):
                f.write(f'{i} {c}\n')

    plot_name = args.plot_name[0] if args.plot_name else 'persistence_length.png'
    pl = fit_PL(correlations, plot_name)
    
    step_nm = l0 * 0.8518
    
    print("\n" + "="*60)
    print("PERSISTENCE LENGTH RESULTS")
    print("="*60)
    print(f"Mean base-pair step length: {step_nm:.4f} nm")
    print(f"\nPersistence length:")
    print(f"  {pl:.1f} base steps")
    print(f"  {pl * step_nm:.1f} nm")
    print("="*60 + "\n")
    print(f"--- {time.time() - start_time:.2f} seconds ---")


if __name__ == '__main__':
    main()

