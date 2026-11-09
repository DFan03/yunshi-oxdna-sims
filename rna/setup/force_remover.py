#!/usr/bin/env python3
import sys
import re

def should_remove(block, ids_to_remove):
    """Check if a block contains particle/ref_particle in ids_to_remove."""
    for line in block:
        match_particle = re.search(r'particle\s*=\s*(\d+)', line)
        match_ref = re.search(r'ref_particle\s*=\s*(\d+)', line)

        if match_particle and int(match_particle.group(1)) in ids_to_remove:
            return True
        if match_ref and int(match_ref.group(1)) in ids_to_remove:
            return True
    return False


def process_file(filename, ids_to_remove):
    with open(filename, 'r') as f:
        lines = f.readlines()

    output = []
    inside_block = False
    current_block = []

    for line in lines:
        stripped = line.strip()

        if stripped.startswith("{"):
            inside_block = True
            current_block = [line]
        elif inside_block:
            current_block.append(line)

            if stripped.endswith("}"):
                # End of block
                if not should_remove(current_block, ids_to_remove):
                    output.extend(current_block)
                # Reset
                inside_block = False
                current_block = []
        else:
            # Outside any block, keep as-is
            output.append(line)

    return output


def main():
    if len(sys.argv) < 3:
        print("Usage: python force_pair_remover.py <force_file> <id1> <id2> ...")
        sys.exit(1)

    filename = sys.argv[1]
    ids_to_remove = set(map(int, sys.argv[2:]))

    cleaned_lines = process_file(filename, ids_to_remove)

    # Print to stdout (so you can redirect if needed)
    for line in cleaned_lines:
        print(line, end='')


if __name__ == "__main__":
    main()