#!/usr/bin/env python3
"""Mean-reference pairwise signed bend and twist analysis.

Usage:
    python pairwise_bend_twist.py traj.dat junctions.txt k average.dat [exclude_fl]

average.dat must be an aligned oxDNA mean structure. The output format is kept
compatible with the signed-bend Jupyter notebook:
    t,bend_signed_deg,twist_deg
"""
import csv
import sys
import numpy as np

EPS = 1e-8
MEAN_BEND_EPS_DEG = 1e-3


def unit(v, label):
    v = np.asarray(v, float)
    n = np.linalg.norm(v)
    if n < EPS:
        raise ValueError(f"Cannot define {label}: nearly zero vector")
    return v / n


def read_junctions(path):
    result = {}
    with open(path, encoding="utf-8") as f:
        for raw in f:
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            p = [x.strip() for x in line.split(",")]
            nums = [int(x) for x in p[1:]]
            ranges = [nums] if len(nums) == 2 else [nums[:2], nums[2:]] if len(nums) == 4 else None
            if ranges is None:
                raise ValueError(f"Expected name,start,end[,start2,end2]: {raw.rstrip()}")
            ids = []
            for start, end in ranges:
                if start < 1 or end < start:
                    raise ValueError(f"Invalid range: {raw.rstrip()}")
                ids.extend(range(start - 1, end))
            result[p[0]] = ids
    return result


def first_frame_size(path):
    started = False
    n = 0
    with open(path, encoding="utf-8") as f:
        for raw in f:
            s = raw.strip()
            if s.startswith("t ="):
                if started:
                    break
                started = True
            elif started and s and not s.startswith("b") and not s.startswith("E"):
                n += 1
    if not n:
        raise ValueError(f"No coordinates in first frame of {path}")
    return n


def frames(path, required):
    required = set(required)
    time, coords, index = None, {}, 0
    with open(path, encoding="utf-8") as f:
        for raw in f:
            s = raw.strip()
            if s.startswith("t ="):
                if time is not None:
                    if required - set(coords):
                        raise ValueError(f"Frame t={time} lacks requested coordinates")
                    yield time, coords
                time, coords, index = float(s.split("=", 1)[1]), {}, 0
            elif not s or s.startswith("b") or s.startswith("E"):
                continue
            else:
                p = s.split()
                if len(p) < 3:
                    raise ValueError(f"Malformed coordinate line: {raw.rstrip()}")
                if index in required:
                    coords[index] = np.array([float(p[0]), float(p[1]), float(p[2])])
                index += 1
    if time is not None:
        if required - set(coords):
            raise ValueError(f"Final frame t={time} lacks requested coordinates")
        yield time, coords


def pts(c, ids): return np.asarray([c[i] for i in ids], float)
def com(c, ids): return pts(c, ids).mean(axis=0)


def markers(c, monomers):
    return np.vstack([pts(c, m[name]) for m in monomers for name in ("LH1", "LH4", "SH1", "SH2")])


def kabsch(ref, cur):
    ref, cur = ref - ref.mean(0), cur - cur.mean(0)
    u, _, vt = np.linalg.svd(ref.T @ cur)
    r = vt.T @ u.T
    if np.linalg.det(r) < 0:
        vt[-1] *= -1
        r = vt.T @ u.T
    return r  # current vector = r @ reference vector


def bend_vector(c, mi, mj):
    ti = unit(com(c, mi["LH4"]) - com(c, mi["LH1"]), "long axis i")
    tj = unit(com(c, mj["LH4"]) - com(c, mj["LH1"]), "long axis j")
    dot = float(np.clip(np.dot(ti, tj), -1, 1))
    sine = float(np.linalg.norm(np.cross(ti, tj)))
    angle = float(np.degrees(np.arctan2(sine, dot)))
    return angle, np.zeros(3) if sine < EPS else angle * unit(tj - dot * ti, "bend direction")


def fallback_direction(avg, mi):
    lh1, lh4 = com(avg, mi["LH1"]), com(avg, mi["LH4"])
    axis = unit(lh4 - lh1, "mean long axis")
    v = com(avg, mi["SH1"]) - 0.5 * (lh1 + lh4)
    return unit(v - np.dot(v, axis) * axis, "mean SH1 direction")


def twist(c, mi, mj):
    lh1, lh4 = com(c, mi["LH1"]), com(c, mj["LH4"])
    sh1, sh2 = com(c, mi["SH1"]), com(c, mj["SH2"])
    axis = unit(lh4 - lh1, "twist axis")
    if np.dot(axis, sh2 - sh1) < 0: axis = -axis
    r1, r2 = sh1 - lh1, sh2 - lh1
    r1 -= np.dot(r1, axis) * axis; r2 -= np.dot(r2, axis) * axis
    r1, r2 = unit(r1, "SH1 radial"), unit(r2, "SH2 radial")
    return float(-np.degrees(np.arctan2(np.dot(axis, np.cross(r1, r2)), np.clip(np.dot(r1, r2), -1, 1))))


def main():
    if len(sys.argv) not in (5, 6):
        print(f"Usage: {sys.argv[0]} traj.dat junctions.txt k average.dat [exclude_fl]", file=sys.stderr); sys.exit(1)
    traj, junction_file, avg_file = sys.argv[1], sys.argv[2], sys.argv[4]
    try:
        k = int(sys.argv[3]); base = read_junctions(junction_file)
        if k < 4: raise ValueError("k must be at least 4")
        missing = [x for x in ("LH1", "LH4", "SH1", "SH2") if x not in base]
        if missing: raise ValueError(f"junctions.txt lacks: {', '.join(missing)}")
        n = first_frame_size(traj)
        if first_frame_size(avg_file) != n or n % k: raise ValueError("trajectory/average size is invalid for k")
        length = n // k
        if max(x for ids in base.values() for x in ids) >= length: raise ValueError("Junction index exceeds monomer length")
    except (OSError, ValueError) as e:
        print(f"Error: {e}", file=sys.stderr); sys.exit(1)

    monomers = [{name: [x + m * length for x in ids] for name, ids in base.items()} for m in range(k)]
    needed = sorted({x for m in monomers for ids in m.values() for x in ids})
    pairs = [(i, i + 1) for i in range(1, k - 2)]
    try:
        _, avg = next(frames(avg_file, needed)); avg_markers = markers(avg, monomers)
        refs = {}
        for i, j in pairs:
            _, v = bend_vector(avg, monomers[i], monomers[j])
            refs[(i, j)] = unit(v, "mean bend direction") if np.linalg.norm(v) >= MEAN_BEND_EPS_DEG else fallback_direction(avg, monomers[i])
    except (ValueError, StopIteration) as e:
        print(f"Error reading average.dat: {e}", file=sys.stderr); sys.exit(1)

    output = {}
    try:
        for i, j in pairs:
            h = open(f"{i+1}{j+1}_pair.csv", "w", newline="", encoding="utf-8")
            w = csv.writer(h); w.writerow(["t", "bend_signed_deg", "twist_deg"]); output[(i, j)] = (h, w)
        for t, c in frames(traj, needed):
            r = kabsch(avg_markers, markers(c, monomers))
            for i, j in pairs:
                try:
                    _, current_v = bend_vector(c, monomers[i], monomers[j])
                    signed = float(np.dot(r.T @ current_v, refs[(i, j)]))
                    tw = twist(c, monomers[i], monomers[j])
                except ValueError:
                    signed, tw = np.nan, np.nan
                output[(i, j)][1].writerow([t, signed, tw])
    finally:
        for h, _ in output.values(): h.close()
    print("Wrote " + ", ".join(f"{i+1}{j+1}_pair.csv" for i, j in pairs))

if __name__ == "__main__":
    main()

