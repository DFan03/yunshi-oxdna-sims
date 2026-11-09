import json
import matplotlib.pyplot as plt

# Load bond instability data
with open("bonds_instability.json", "r") as f:
    data = json.load(f)

instability = data["bond_instability_index"]

# Generate X-axis (nucleotide index)
x = list(range(1, len(instability) + 1))

# Create the plot
plt.figure(figsize=(12, 4))  # wide figure for sequences
plt.plot(x, instability, color='red', linewidth=1.5)
plt.xlabel("Nucleotide index", fontsize=12)
plt.ylabel("Bond Instability Index", fontsize=12)
plt.title("Per-nucleotide Bond Instability", fontsize=14)
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()

# Save as JPEG with 300 dpi
plt.savefig("bond_instability_plot.jpeg", dpi=300, format='jpeg')
plt.close()

print("Saved bond instability plot as bond_instability_plot.jpeg at 300 dpi")
