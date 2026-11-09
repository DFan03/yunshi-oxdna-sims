#!/usr/bin/env python3
import os
import shutil
import glob
import re

# ---- Step 1: Copy files ----

src_input = "../input_pre_relax"
src_conf = "../last_conf.dat"

if os.path.exists(src_input):
    shutil.copy(src_input, ".")
else:
    print(f"Warning: {src_input} not found")

if os.path.exists(src_conf):
    shutil.copy(src_conf, ".")
else:
    print(f"Warning: {src_conf} not found")

# copy all .top files from parent
top_files = glob.glob("../*.top")
for f in top_files:
    shutil.copy(f, ".")

if not top_files:
    print("Warning: No .top files found in parent directory")

# use the first .top file for replacement
top_name = os.path.basename(top_files[0]) if top_files else None


# ---- Step 2: Modify *input.txt ----

input_files = glob.glob("*input.txt")

for infile in input_files:
    with open(infile, "r") as f:
        text = f.read()

    if top_name:
        text = re.sub(r"output\.top", top_name, text)

    text = re.sub(r"output\.dat", "last_conf.dat", text)

    with open(infile, "w") as f:
        f.write(text)

    print(f"Updated {infile}")


# ---- Step 3: Generate cmds file ----

cmds = [
    "oat generate_forces input_pre_relax last_conf.dat",
    "oat forces2pairs forces.txt",
    "oat bond_analysis *input.txt *_trj.dat pairs.txt -p 4",
    "oat mean -p 4 -o mean.dat -d devs.json *_trj.dat"
]

with open("cmds.sh", "w") as f:
    f.write("#!/bin/bash\n\n")
    for c in cmds:
        f.write(c + "\n")

os.chmod("cmds.sh", 0o755)

print("Created cmds.sh")
