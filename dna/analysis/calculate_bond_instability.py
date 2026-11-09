import json

input_file = "bonds.json"
output_file = "bonds_instability.json"

# Load data
with open(input_file, "r") as f:
    data = json.load(f)

occupancy = data["occupancy"]

# Transform to Bond Instability Index
instability = [4 * p * (1 - p) for p in occupancy]

# Save result
out = {"bond_instability_index": instability}

with open(output_file, "w") as f:
    json.dump(out, f, indent=2)

print(f"Saved transformed data to {output_file}")
