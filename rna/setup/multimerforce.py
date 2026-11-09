#!/usr/bin/env python3

import argparse

MONOMER_LENGTH = 244

# Local positions within one monomer (0–243)
A_POS = [91, 90, 89, 88, 87, 86]
B_POS = [152, 153, 154, 155, 156, 157]
C_POS = [58, 57, 56, 55, 54, 53]
D_POS = [185, 186, 187, 188, 189, 190]


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
    Generate all base-pair index pairs.

    Linear (circular = False):
      - a_0 and c_0 are unpaired
      - b_{N-1} and d_{N-1} are unpaired
      - for i = 1..N-1:
          a_i <-> b_{i-1}
          c_i <-> d_{i-1}

    Circular (circular = True):
      - for i = 0..N-1, with prev = (i-1) mod N:
          a_i <-> b_prev
          c_i <-> d_prev
    """
    pairs = []

    if n_oligos <= 1:
        return pairs

    if circular:
        # Full ring: connect every chain to its previous neighbor
        for i in range(n_oligos):
            prev = (i - 1) % n_oligos
            offset_i = i * MONOMER_LENGTH
            offset_prev = prev * MONOMER_LENGTH

            # a_i with b_prev
            for a_local, b_local in zip(A_POS, B_POS):
                pairs.append((offset_i + a_local, offset_prev + b_local))

            # c_i with d_prev
            for c_local, d_local in zip(C_POS, D_POS):
                pairs.append((offset_i + c_local, offset_prev + d_local))
    else:
        # Linear: skip i = 0 so a_0 and c_0 are unpaired; last b,d remain unpaired
        for i in range(1, n_oligos):
            prev = i - 1
            offset_i = i * MONOMER_LENGTH
            offset_prev = prev * MONOMER_LENGTH

            # a_i with b_{i-1}
            for a_local, b_local in zip(A_POS, B_POS):
                pairs.append((offset_i + a_local, offset_prev + b_local))

            # c_i with d_{i-1}
            for c_local, d_local in zip(C_POS, D_POS):
                pairs.append((offset_i + c_local, offset_prev + d_local))

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
            "and write corresponding mutual_trap forces to forces.txt."
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
