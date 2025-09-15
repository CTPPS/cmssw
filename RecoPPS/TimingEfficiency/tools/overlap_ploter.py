import json
import matplotlib.pyplot as plt
import os
import argparse

# -------------------------
# Argument parser
# -------------------------
parser = argparse.ArgumentParser(
    description="Plot merged JSON data for specific steps and types."
)
parser.add_argument(
    "--n", type=int, default=1, help="Step size for x-ticks (default: 1)."
)
parser.add_argument(
    "--json1", type=str, required=True, help="Path to the first JSON file (all bunches)."
)
parser.add_argument(
    "--json2", type=str, required=True, help="Path to the second JSON file (picked bunches)."
)
parser.add_argument(
    "--range",
    type=int,
    nargs=2,
    metavar=("START", "END"),
    help="Range of runs (inclusive) to plot.",
)
parser.add_argument(
    "--width", type=float, default=15, help="Figure width (default: 15)."
)
parser.add_argument(
    "--output",
    type=str,
    default="plots",
    help="Output directory for plots (default: plots).",
)

args = parser.parse_args()

# -------------------------
# Load JSON data
# -------------------------
with open(args.json1, "r") as f:
    data = json.load(f)

with open(args.json2, "r") as f:
    data2 = json.load(f)


# Create output folder for plots
output_dir = args.output
os.makedirs(output_dir, exist_ok=True)

# Sort and filter keys
all_keys = sorted(data.keys())
if args.range:
    start, end = args.range
    sorted_keys = [k for k in all_keys if start <= int(k) <= end]
else:
    sorted_keys = all_keys

x = [int(k) for k in sorted_keys]
x_positions = range(len(sorted_keys))

num_of_diams = {
    "box": max(
        max(len(data[k]["45"]["box"]), len(data[k]["56"]["box"])) for k in sorted_keys
    ),
    "cyl": max(
        max(len(data[k]["45"]["cyl"]), len(data[k]["56"]["cyl"])) for k in sorted_keys
    ),
}


# -------------------------
# Extract values function
# -------------------------
def extract_values(step, type_key):
    return [
        [data[k][step][type_key][i] for k in sorted_keys]
        for i in range(num_of_diams[type_key])
    ]


def extract_values2(step, type_key):
    return [
        [data2[k][step][type_key][i] for k in sorted_keys]
        for i in range(num_of_diams[type_key])
    ]


# -------------------------
# Plot and save
# -------------------------
def plot_group(step, type_key):
    values = extract_values(step, type_key)
    values2 = extract_values2(step, type_key)
    fig, axs = plt.subplots(
        num_of_diams[type_key],
        1,
        figsize=(args.width, 3.5 * num_of_diams[type_key]),
        sharex=True,
    )
    fig.suptitle(f"{step} {type_key}")

    for i in range(num_of_diams[type_key]):
        axs[i].scatter(x_positions, values[i], marker="o")
        axs[i].set_xticks(x_positions[:: args.n])
        axs[i].set_xticklabels(sorted_keys[:: args.n], rotation=45, ha="right")
        axs[i].set_ylabel(f"Diam {i+1} efficiency")
        axs[i].set_ylim(0, 1.05)
        axs[i].grid(True)

        axs[i].scatter(
            x_positions,
            values2[i],
            marker="o",
            facecolors="none",  # empty
            edgecolors="red",  # red border
            linewidth=1.5,
            label="Picked bunches",
        )

    axs[num_of_diams[type_key] - 1].set_xlabel("Run", labelpad=15)
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    # Save the figure as PNG
    filename = os.path.join(output_dir, f"{step}_{type_key}.png")
    plt.savefig(filename)
    plt.close(fig)


# -------------------------
# Generate and save plots
# -------------------------
plot_group("45", "box")
plot_group("45", "cyl")
plot_group("56", "box")
plot_group("56", "cyl")

print(f"Plots saved in folder: {output_dir}")
