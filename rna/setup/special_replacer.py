import glob

# Find the .oxdna and .top files
oxdna_files = glob.glob("*.dat")
top_files = glob.glob("*.top")

if not oxdna_files or not top_files:
    raise FileNotFoundError("Missing .oxdna or .top file in directory")

oxdna_name = oxdna_files[0]
top_name = top_files[0]

# Process all input* files
for fname in glob.glob("input*"):
    with open(fname, "r") as f:
        lines = f.readlines()

    new_lines = []
    for line in lines:
        line = line.replace("replace_me.dat", oxdna_name)
        line = line.replace("replace_me.top", top_name)
        new_lines.append(line)

    with open(fname, "w") as f:
        f.writelines(new_lines)

print("Replacement completed.")
