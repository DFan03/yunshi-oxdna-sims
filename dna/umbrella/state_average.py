#!/usr/bin/env python3

"""
Compute weighted average free energy (kBT) for open and closed states.

Inputs (in current directory):

  open_closed_profile.dat :
    # x_sim  x_nm  P_r  P_open_given_r  P_closed_given_r  N_tot  N_open  N_closed
    ...

  com_wham_profile.dat :
    # x_sim  x_nm  F(kBT)  P
    ...

For each row in open_closed_profile.dat:

  1. Try to find a matching WHAM row by x_nm with decreasing decimal precision:
       - First try highest precision (up to 6 decimals, e.g. 6.207055).
       - Then 5 decimals (6.20706), 4 decimals (6.2071), ...
       - Down to 1 decimal (6.2).
     Matching rule at precision d:
       round(x_nm_open_closed, d) == round(x_nm_wham, d)
     If multiple WHAM rows match at that precision, the one with
     smallest absolute |x_nm_open_closed - x_nm_wham| is chosen.

  2. For each row:
       - Print: distance (x_nm), whether a match was found,
         the decimal precision used, and the matched WHAM x_nm.

  3. If no WHAM row matches down to 1 decimal:
       - Print a warning for that row.
       - Ignore that row in the averages.

  4. For matched rows:
       - Let F be F(kBT) from WHAM.
       - Accumulate:
           kbt_open_total   += F * N_open
           kbt_closed_total += F * N_close
           N_open_total     += N_open
           N_closed_total   += N_close

Finally:

  F_open_avg   = kbt_open_total   / N_open_total   (if N_open_total   > 0)
  F_closed_avg = kbt_closed_total / N_closed_total (if N_closed_total > 0)
"""

OPEN_CLOSED_FILE = "open_closed_profile.dat"
WHAM_FILE        = "com_wham_profile.dat"

# Maximum number of decimal places we try for matching
MAX_DECIMALS = 6
MIN_DECIMALS = 1


def load_wham(wham_file):
    """
    Load WHAM data.

    Returns:
      wham_rows : list of (x_nm, F_kBT)
    """
    wham_rows = []
    with open(wham_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            # columns: x_sim, x_nm, F(kBT), P
            x_nm  = float(parts[1])
            F_kBT = float(parts[2])
            wham_rows.append((x_nm, F_kBT))
    if not wham_rows:
        raise RuntimeError(f"No data rows found in WHAM file: {wham_file}")
    return wham_rows


def find_best_match(x_nm_target, wham_rows, max_decimals=MAX_DECIMALS, min_decimals=MIN_DECIMALS):
    """
    Try to match x_nm_target to WHAM x_nm with decreasing decimal precision.

    At each precision d (from max_decimals down to min_decimals):
      - Round target: rt = round(x_nm_target, d)
      - Find all WHAM rows whose round(x_nm_wham, d) == rt
      - If one or more candidates exist:
          * pick the one with smallest |x_nm_target - x_nm_wham|
          * return (F_kBT, used_decimals, matched_x_nm)

    If no match found down to min_decimals:
      return (None, 0, None)
    """
    best_match = None
    best_diff  = None
    best_dec   = 0

    for dec in range(max_decimals, min_decimals - 1, -1):
        rt = round(x_nm_target, dec)

        candidates = []
        for x_nm_wham, F_kBT in wham_rows:
            if round(x_nm_wham, dec) == rt:
                candidates.append((x_nm_wham, F_kBT))

        if not candidates:
            continue

        # Choose candidate with smallest absolute difference
        x_nm_best, F_best = min(
            candidates,
            key=lambda row: abs(row[0] - x_nm_target)
        )
        diff = abs(x_nm_best - x_nm_target)

        best_match = (x_nm_best, F_best)
        best_diff  = diff
        best_dec   = dec
        # As soon as we find a match at this precision, we can break;
        # we don't try lower precisions.
        break

    if best_match is None:
        return None, 0, None

    return best_match[1], best_dec, best_match[0]  # F_kBT, dec, matched_x_nm


def compute_weighted_averages(open_closed_file, wham_file):
    wham_rows = load_wham(wham_file)

    kbt_open_total   = 0.0
    kbt_closed_total = 0.0
    N_open_total     = 0
    N_closed_total   = 0

    with open(open_closed_file) as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            # x_sim  x_nm  P_r  P_open_given_r  P_closed_given_r  N_tot  N_open  N_closed
            if len(parts) < 8:
                continue

            x_nm_str = parts[1]
            x_nm     = float(x_nm_str)
            N_open   = int(parts[6])
            N_close  = int(parts[7])

            # Skip bins with no open/closed counts
            if N_open == 0 and N_close == 0:
                print(
                    f"[row {lineno}] x_nm={x_nm_str} -> "
                    f"no open/closed counts (N_open=N_close=0); ignored."
                )
                continue

            F_kBT, used_dec, matched_x_nm = find_best_match(x_nm, wham_rows)

            if F_kBT is None:
                print(
                    f"[row {lineno}] x_nm={x_nm_str} -> "
                    f"no WHAM match down to {MIN_DECIMALS} decimal(s); row ignored."
                )
                continue

            print(
                f"[row {lineno}] x_nm={x_nm_str} -> MATCH: "
                f"WHAM x_nm={matched_x_nm:.6f}, decimals={used_dec}"
            )

            # Accumulate weighted sums
            kbt_open_total   += F_kBT * N_open
            kbt_closed_total += F_kBT * N_close
            N_open_total     += N_open
            N_closed_total   += N_close

    F_open_avg   = kbt_open_total   / N_open_total   if N_open_total   > 0 else None
    F_closed_avg = kbt_closed_total / N_closed_total if N_closed_total > 0 else None

    return F_open_avg, F_closed_avg, N_open_total, N_closed_total


def main():
    F_open_avg, F_closed_avg, N_open_total, N_closed_total = compute_weighted_averages(
        OPEN_CLOSED_FILE, WHAM_FILE
    )

    print("\n=== Weighted average free energies (k_B T units) ===")
    print(f"N_open_total   = {N_open_total}")
    print(f"N_closed_total = {N_closed_total}")

    if F_open_avg is not None:
        print(f"F_open_avg     = {F_open_avg:.6f} k_B T")
    else:
        print("F_open_avg     = N/A (no open frames)")

    if F_closed_avg is not None:
        print(f"F_closed_avg   = {F_closed_avg:.6f} k_B T")
    else:
        print("F_closed_avg   = N/A (no closed frames)")


if __name__ == "__main__":
    main()
