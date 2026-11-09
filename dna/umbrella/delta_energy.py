#!/usr/bin/env python3

"""
Compute KL-open / KL-closed state probabilities and ΔG from WHAM + HB fractions.

This implements exactly what was requested:

  - Define KL-open  : HB = 0
  - Define KL-closed: HB = 6

  - For each bin (window):
      * Take its unbiased probability P_r from WHAM (PMF/reweighting).
      * Multiply by the fraction of frames in that bin in each HB class:
            P(open | r)   and  P(closed | r)
      * Sum these contributions over all bins:

          P_open_tilde   = Σ_j P_r(j) * P(open | r_j)
          P_closed_tilde = Σ_j P_r(j) * P(closed | r_j)

  - Restrict to the {open, closed} subspace:
        P_open   = P_open_tilde   / (P_open_tilde + P_closed_tilde)
        P_closed = P_closed_tilde / (P_open_tilde + P_closed_tilde)

  - Free-energy difference:
        ΔG_open→closed = -k_B T ln(P_closed / P_open)

Input format (open_closed_profile.dat):

  # x_sim  x_nm  P_r  P_open_given_r  P_closed_given_r  N_tot  N_open  N_closed
    ...

Here:
  - P_r              : unbiased P(r) from WHAM
  - P_open_given_r   : fraction of open (HB=0) frames in that bin
  - P_closed_given_r : fraction of closed (HB=6) frames in that bin
"""

import math

OPEN_CLOSED_FILE = "open_closed_profile.dat"

# Temperature + gas constant (for kJ/mol output, optional)
TEMPERATURE_K = 293.0
R_KJMOLK      = 8.314e-3  # kJ/(mol·K)
RT_KJMOL      = R_KJMOLK * TEMPERATURE_K


def compute_state_probabilities(open_closed_file):
    """
    Read open_closed_profile.dat and compute:

      - P_open, P_closed (normalized over {open, closed})
      - ΔG_open→closed in k_B T and kJ/mol

    Uses:
      P_open_tilde   = Σ P_r * P_open_given_r
      P_closed_tilde = Σ P_r * P_closed_given_r

    Bins with N_tot == 0 are ignored.
    """
    P_open_tilde   = 0.0
    P_closed_tilde = 0.0

    total_bins = 0
    used_bins  = 0

    with open(open_closed_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            total_bins += 1

            parts = line.split()
            # x_sim  x_nm  P_r  P_open_given_r  P_closed_given_r  N_tot  N_open  N_closed
            if len(parts) < 8:
                continue

            P_r            = float(parts[2])
            P_open_given_r = float(parts[3])
            P_closed_given_r = float(parts[4])
            N_tot          = int(parts[5])

            # Skip bins with no statistics
            if N_tot <= 0:
                continue

            used_bins += 1

            P_open_tilde   += P_r * P_open_given_r
            P_closed_tilde += P_r * P_closed_given_r

    if used_bins == 0:
        raise RuntimeError("No bins with N_tot > 0 found in open_closed_profile.dat")

    P_sum = P_open_tilde + P_closed_tilde
    if P_sum <= 0.0:
        raise RuntimeError(
            "Total reweighted probability for open+closed is zero; "
            "check input data or HB definitions."
        )

    P_open   = P_open_tilde   / P_sum
    P_closed = P_closed_tilde / P_sum

    # ΔG in k_B T
    dG_kBT = -math.log(P_closed / P_open)

    # ΔG in kJ/mol (optional)
    dG_kJmol = dG_kBT * RT_KJMOL

    return {
        "P_open": P_open,
        "P_closed": P_closed,
        "dG_kBT": dG_kBT,
        "dG_kJmol": dG_kJmol,
        "P_open_tilde": P_open_tilde,
        "P_closed_tilde": P_closed_tilde,
        "total_bins": total_bins,
        "used_bins": used_bins,
    }


def main():
    res = compute_state_probabilities(OPEN_CLOSED_FILE)

    print("=== KL-open / KL-closed state probabilities (from WHAM + HB fractions) ===")
    print(f"Total bins read    : {res['total_bins']}")
    print(f"Bins used (N_tot>0): {res['used_bins']}")
    print("")
    print(f"P_open_tilde       = {res['P_open_tilde']:.6e}")
    print(f"P_closed_tilde     = {res['P_closed_tilde']:.6e}")
    print(f"P_open (norm.)     = {res['P_open']:.6f}")
    print(f"P_closed (norm.)   = {res['P_closed']:.6f}")
    print("")
    print(f"ΔG_open→closed     = {res['dG_kBT']:.6f} k_B T")
    print(f"                   ≈ {res['dG_kJmol']:.6f} kJ/mol (T = {TEMPERATURE_K} K)")


if __name__ == "__main__":
    main()
