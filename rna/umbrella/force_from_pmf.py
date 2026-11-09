#!/usr/bin/env python3
"""
compute_force_from_wham.py

Reads:
  com_wham_profile.dat  (produced by wham_com_completed.py)
and computes:
  - f(x) ≈ -dF/dx   as a function of x (nm), in pN
  - average force between first and last window using ΔF/Δx

Assumptions:
  - Column 1: x_nm   (nm)
  - Column 2: F      (k_B T)
"""

import numpy as np

# Conversion: 1 k_B T / nm ≈ 4.1 pN (at room T ~293 K)
KBT_PER_NM_TO_PN = 4.1

INFILE  = "com_wham_profile.dat"
OUTFILE = "com_force_profile.dat"

def main():
    print("=== compute_force_from_wham.py ===")
    print(f"Reading WHAM profile from {INFILE}")

    # Load data, skipping comment lines
    data = np.loadtxt(INFILE, comments="#")
    # According to wham_com_completed.py:
    # col0: x_sim, col1: x_nm, col2: F, col3: P
    x_nm = data[:, 1]
    F    = data[:, 2]

    # Numerical derivative dF/dx with central differences
    # f(x) ≈ -dF/dx * 4.1 pN
    dFdx = np.zeros_like(F)

    # central differences for interior points
    dFdx[1:-1] = (F[2:] - F[:-2]) / (x_nm[2:] - x_nm[:-2])

    # forward/backward for endpoints
    dFdx[0] = (F[1] - F[0]) / (x_nm[1] - x_nm[0])
    dFdx[-1] = (F[-1] - F[-2]) / (x_nm[-1] - x_nm[-2])

    force_pN = -dFdx * KBT_PER_NM_TO_PN

    # Save force profile
    print(f"Writing force–extension profile to {OUTFILE}")
    with open(OUTFILE, "w") as f:
        f.write("# x_nm  F_kBT  dFdx_kBT_per_nm  force_pN\n")
        for xi, Fi, dFi, fi in zip(x_nm, F, dFdx, force_pN):
            f.write(f"{xi:.6f}  {Fi:.6f}  {dFi:.6f}  {fi:.6f}\n")

    # Simple average force over full span (first -> last)
    dF_total = F[-1] - F[0]
    dx_total = x_nm[-1] - x_nm[0]
    avg_force_pN = - (dF_total / dx_total) * KBT_PER_NM_TO_PN

    print("\n=== Summary over full profile range ===")
    print(f"x_first ≈ {x_nm[0]:.3f} nm, F_first ≈ {F[0]:.3f} k_B T")
    print(f"x_last  ≈ {x_nm[-1]:.3f} nm, F_last  ≈ {F[-1]:.3f} k_B T")
    print(f"ΔF = F_last - F_first = {dF_total:.3f} k_B T")
    print(f"Δx = x_last - x_first = {dx_total:.3f} nm")
    print(f"Average force ≈ {avg_force_pN:.3f} pN (using -ΔF/Δx)")
    print("=======================================")


if __name__ == "__main__":
    main()
