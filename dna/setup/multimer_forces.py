#!/usr/bin/env python3

import argparse

MONOMER_LENGTH = 258

# Local positions within one monomer (0–243)
A_POS = [109,108,107,106,105,104]
B_POS = [170,171,172,173,174,175]
C_POS = [208,207,206,205,204,203]
D_POS = [71,72,73,74,75,76]


def parse_bool(value: str) -> bool:
    v = value.strip().lower()
    if v in ("true", "t", "1", "yes", "y"):
        return True
    if v in ("false", "f", "0", "no", "n"):
        return False
    raise argparse.ArgumentTypeError(
        f"Cannot interpret '{value}' as boolean; use true/false."
    )


def generate_pairs(n_oligos: int, circular: bool):
    """
    Generate all base-pair index pairs for the NEW geometry.

    Linear (circular = False):
      - A_0 and D_0 are unpaired
      - B_{N-1} and C_{N-1} are unpaired
      - for i = 0..N-2:
          B_i <-> A_{i+1}
          C_i <-> D_{i+1}

    Circular (circular = True) [assumed]:
      - for i = 0..N-1, with nxt = (i+1) mod N:
          B_i <-> A_nxt
          C_i <-> D_nxt
    """
    pairs = []

    if n_oligos <= 1:
        # No inter-monomer contacts in both cases
        return pairs

    if circular:
        # Full ring: connect every chain to its next neighbor
        for i in range(n_oligos):
            nxt = (i + 1) % n_oligos
            offset_i = i * MONOMER_LENGTH
            offset_nxt = nxt * MONOMER_LENGTH

            # B_i with A_{nxt}
            for b_local, a_local in zip(B_POS, A_POS):
                pairs.append((offset_i + b_local, offset_nxt + a_local))

            # C_i with D_{nxt}
            for c_local, d_local in zip(C_POS, D_POS):
                pairs.append((offset_i + c_local, offset_nxt + d_local))
    else:
        # Linear: i goes up to N-2; last monomer's B,C unpaired
        for i in range(n_oligos - 1):
            nxt = i + 1
            offset_i = i * MONOMER_LENGTH
            offset_nxt = nxt * MONOMER_LENGTH

            # B_i with A_{i+1}
            for b_local, a_local in zip(B_POS, A_POS):
                pairs.append((offset_i + b_local, offset_nxt + a_local))

            # C_i with D_{i+1}
            for c_local, d_local in zip(C_POS, D_POS):
                pairs.append((offset_i + c_local, offset_nxt + d_local))

    return pairs


def write_forces_file(pairs, stiffness: float, filename: str = "forces.txt"):
    """
    Write mutual_trap force blocks for each pair into forces.txt.

    For each (i, j) in pairs, write TWO blocks:
      1) particle = i, ref_particle = j
      2) particle = j, ref_particle = i
    """
    with open(filename, "w") as f:
        for i, j in pairs:
            # i -> j
            f.write("{\n")
            f.write("type = mutual_trap\n")
            f.write(f"particle = {i}\n")
            f.write(f"stiff = {stiffness}\n")
            f.write("r0 = 1.2\n")
            f.write(f"ref_particle = {j}\n")
            f.write("PBC=1\n")
            f.write("}\n")

            # j -> i
            f.write("{\n")
            f.write("type = mutual_trap\n")
            f.write(f"particle = {j}\n")
            f.write(f"stiff = {stiffness}\n")
            f.write("r0 = 1.2\n")
            f.write(f"ref_particle = {i}\n")
            f.write("PBC=1\n")
            f.write("}\n")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Generate base-pair indices for an RNA origami homo-oligomer\n"
            "with B,C -> next A,D geometry and write mutual_trap forces."
        )
    )
    parser.add_argument(
        "n_oligos",
        type=int,
        help="Number of oligomer units (monomers).",
    )
    parser.add_argument(
        "circular",
        type=parse_bool,
        help="Whether the oligomer is circular (true/false).",
    )
    parser.add_argument(
        "stiffness",
        type=float,
        help="Stiffness value to use in the mutual_trap forces.",
    )

    args = parser.parse_args(argv)

    if args.n_oligos < 1:
        parser.error("n_oligos must be >= 1")

    pairs = generate_pairs(args.n_oligos, args.circular)
    write_forces_file(pairs, args.stiffness)


if __name__ == "__main__":
    main()
