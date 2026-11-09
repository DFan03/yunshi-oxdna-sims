#!/usr/bin/env python3
"""
compute_barrier_forces.py

Reads:
  com_wham_profile.dat  (produced by wham_com_completed.py)

Assumptions:
  - COM distance x_nm on the RIGHT (larger) corresponds to OPEN state.
  - COM distance x_nm on the LEFT (smaller) corresponds to CLOSED state.

Finds:
  - Barrier as the global maximum of F(x)
  - CLOSED minimum: global minimum on the left of the barrier
  - OPEN   minimum: global minimum on the right of the barrier

Computes:
  - Barrier heights ΔF_closed^‡, ΔF_open^‡ in k_B T
  - Distances Δx_closed^‡, Δx_open^‡ in nm
  - Average downhill forces from barrier to each minimum, in pN

Interpretation:
  - These are average equilibrium driving forces along the COM coordinate,
    after crossing the barrier, toward the CLOSED and OPEN minima.
"""

import numpy as np

KBT_PER_NM_TO_PN = 4.1
INFILE = "com_wham_profile.dat"

def main():
    print("=== compute_barrier_forces.py ===")
    print(f"Reading WHAM profile from {INFILE}")

    data = np.loadtxt(INFILE, comments="#")
    # col0: x_sim, col1: x_nm, col2: F, col3: P
    x_nm = data[:, 1]
    F    = data[:, 2]

    if len(F) < 3:
        print("ERROR: Not enough points in PMF to define minima and barrier.")
        return

    # 1) Barrier = global maximum of F(x)
    i_barrier = int(np.argmax(F))
    x_barrier = x_nm[i_barrier]
    F_barrier = F[i_barrier]

    # Need at least one point on each side to define closed/open minima
    if i_barrier == 0 or i_barrier == len(F) - 1:
        print("WARNING: Global maximum is at the boundary of the PMF.")
        print("         This suggests there is no well-defined interior barrier")
        print("         between two minima in the sampled range.")
        print("\nPMF summary:")
        print(f"  x_min  = {x_nm[0]:.3f} nm, F_min = {F[0]:.3f} k_B T")
        print(f"  x_max  = {x_nm[-1]:.3f} nm, F_max = {F[-1]:.3f} k_B T")
        print(f"  barrier at edge: x = {x_barrier:.3f} nm, F = {F_barrier:.3f} k_B T (index {i_barrier})")
        return

    # 2) CLOSED min = lowest F(x) on the LEFT of the barrier (smaller x)
    F_left  = F[:i_barrier]
    x_left  = x_nm[:i_barrier]
    if len(F_left) == 0:
        print("ERROR: No points to the left of the barrier; cannot define CLOSED minimum.")
        return
    i_min_left_local = int(np.argmin(F_left))
    i_closed = i_min_left_local  # global index

    # 3) OPEN min = lowest F(x) on the RIGHT of the barrier (larger x)
    F_right = F[i_barrier+1:]
    x_right = x_nm[i_barrier+1:]
    if len(F_right) == 0:
        print("ERROR: No points to the right of the barrier; cannot define OPEN minimum.")
        return
    i_min_right_local = int(np.argmin(F_right))
    i_open = i_barrier + 1 + i_min_right_local

    x_closed  = x_nm[i_closed]
    x_open    = x_nm[i_open]
    F_closed  = F[i_closed]
    F_open    = F[i_open]

    # 4) Barrier heights and distances (from barrier to each minimum)
    dF_closed = F_barrier - F_closed
    dF_open   = F_barrier - F_open
    dx_closed = abs(x_barrier - x_closed)
    dx_open   = abs(x_barrier - x_open)

    f_closed_pN = (dF_closed / dx_closed) * KBT_PER_NM_TO_PN if dx_closed > 0 else 0.0
    f_open_pN   = (dF_open   / dx_open)   * KBT_PER_NM_TO_PN if dx_open   > 0 else 0.0

    print("\n=== Identified states (global max + left/right minima) ===")
    print(f"CLOSED minimum (left of barrier): x ≈ {x_closed:.3f} nm,  F ≈ {F_closed:.3f} k_B T (index {i_closed})")
    print(f"OPEN   minimum (right of barrier): x ≈ {x_open:.3f} nm,    F ≈ {F_open:.3f} k_B T (index {i_open})")
    print(f"Barrier (global max):              x ≈ {x_barrier:.3f} nm, F ≈ {F_barrier:.3f} k_B T (index {i_barrier})")

    print("\n=== Barrier properties (downhill from barrier) ===")
    print("From barrier to CLOSED (left):")
    print(f"  ΔF_closed^‡ = F_barrier - F_closed = {dF_closed:.3f} k_B T")
    print(f"  Δx_closed^‡ = |x_barrier - x_closed| = {dx_closed:.3f} nm")
    print(f"  Average downhill closing force ≈ {f_closed_pN:.3f} pN (toward CLOSED)")

    print("\nFrom barrier to OPEN (right):")
    print(f"  ΔF_open^‡ = F_barrier - F_open = {dF_open:.3f} k_B T")
    print(f"  Δx_open^‡ = |x_barrier - x_open| = {dx_open:.3f} nm")
    print(f"  Average downhill opening force ≈ {f_open_pN:.3f} pN (toward OPEN)")

    print("\nInterpretation:")
    print("  These are average equilibrium driving forces along the COM coordinate,")
    print("  after crossing the barrier, pushing the system down toward the CLOSED")
    print("  and OPEN minima. You can call them 'average downhill closing force'")
    print("  and 'average downhill opening force' for your nanotweezer.\n")


if __name__ == "__main__":
    main()
