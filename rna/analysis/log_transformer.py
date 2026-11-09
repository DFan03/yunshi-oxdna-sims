import json
import numpy as np

# Input and output file names
input_file = "bonds_instability.json"
output_file = "log_transformed_bond_instability.json"

# Load the bond instability data
with open(input_file, "r") as f:
    data = json.load(f)

BII = np.array(data["bond_instability_index"])

# Small epsilon to avoid log(0)
epsilon = 1e-6

# Log transform
B_log = np.log(BII + epsilon)

# Min-max scaling to 0-1
B_scaled = (B_log - B_log.min()) / (B_log.max() - B_log.min())

# Save the transformed data
out = {"log_transformed_bond_instability": B_scaled.tolist()}
with open(output_file, "w") as f:
    json.dump(out, f, indent=2)

print(f"Saved log-transformed bond instability to {output_file}")
