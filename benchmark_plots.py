import sys
import numpy as np
import matplotlib.pyplot as plt


semi_rings = {
    "min_plus": 0,
    "or_and":   1,
    "max_min":  2
}

implementations = {
    "baseline":                                 0,
    "basic_opt":                                1,
    "tiled_floyd.c":                            2,
    "vectorized_tiled_floyd.c":                 3,
    "min_plus_generated_vectorized_tiled.c":    4,
    "max_min_generated_vectorized_tiled.c":     4,
    "or_and_generated_vectorized_tiled.c":      4
}

x = [8, 24, 72, 216, 648, 1944]  # N values
y = np.zeros((len(semi_rings), 5, len(x)))     # data of baseline, basic opt, tiled, tiled vectorized
cur_data = []
impl = 0
fwi = 0

if len(sys.argv) != 2:
    print("Error: Results file needed")

# Read plot data from csv file
config_file: str = sys.argv[1]
with open(config_file, "r") as config:

    lines = config.readlines()
    for line in lines:
        # Handle each semi-ring for each implementation
        content = line.split(",")
        if len(content) < len(semi_rings):
            # Skip first pass
            if len(cur_data) > 0:
                y[semi_rings[fwi]][implementations[impl]] = cur_data
                cur_data = []

            (fwi, impl) = content[0].rstrip("\n").split(" ")
        else:
            (_, _, perf) = content
            cur_data.append(perf)

    # Handle last line
    y[semi_rings[fwi]][implementations[impl]] = cur_data


fw = ["(min, +)", "(or, and)", "(max, min)"]
fw_file_name = ["min_plus", "or_and", "max_min"]
vect_peak_performance = 9.94
scalar_peak_performance = 2.48

# Plot the results
for fwi in range(3):
    plt.plot(x, y[fwi][0], label="baseline", marker='o')
    # plt.plot(x, y[fwi][1], label="basic optimizations", marker='o')
    plt.plot(x, y[fwi][2], label="tiled", marker='o')
    plt.plot(x, y[fwi][3], label="tiled vectorized", marker='o')
    plt.plot(x, y[fwi][4], label="generated tiled vectorized", marker='o')

    # plt.axhline(y=vect_peak_performance, label="vectorized peak performance", color='r', linestyle='--')
    plt.axhline(y=scalar_peak_performance, label="scalar peak performance", color='y', linestyle='--')

    # Gray background lines
    for i in range(0, 11):
        plt.axhspan(i, i+.5, facecolor='0.1', alpha=0.1)

    plt.suptitle(fw[fwi] + " performances", fontweight="bold", fontsize=12, x=0.31)
    plt.title("flags: -march=native -O3 -ffast-math", loc='left')
    plt.legend()
    plt.xscale('log', base=3)
    plt.ylim([0.0, 5.5])
    plt.xlabel("N")
    if fwi == 1:
        plt.ylabel("ops/cycle")
    else:
        plt.ylabel("flops/cycle")

    plt.savefig("./plots/" + fw_file_name[fwi] + "_performance_plot.pdf", bbox_inches="tight")
    plt.cla()
