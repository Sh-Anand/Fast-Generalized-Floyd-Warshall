import os
import sys
from subprocess import run

THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
fw = ["min_plus", "max_min", "or_and"]  # functions in generalized floyd warshall
perf_command = "sudo perf stat -B -e cache-references,cache-misses,cycles,instructions,branches,faults,migrations"
flags = " ../tsc_x86.h -march=native -O3 -ffast-math"


def read_config_file():
    # Check for config file
    if len(sys.argv) != 2:
        print("Error: Configuration file needed")
        return

    fwi = -1
    config_file: str = sys.argv[1]
    with open(config_file, "r") as config:
        lines = config.readlines()
        for line in lines:
            # Handle each semi-ring
            if len(line.split(",")) < 2:
                fwi = fwi + 1
                continue

            (n, l1, b, _) = line.split(",")
            print("Begin benchmarking " + str(fw[fwi]) + " with N = %s, L1 = %s, and B = %s" % (n, l1, b))

            if int(n) < 216:     # If small N, run without unrolling
                file_abs_path = os.path.join(THIS_FOLDER, "../vectorized_tiled_floyd.c")
                run("gcc -o ffw " + file_abs_path + flags, shell=True)
                executable_abs_path = os.path.join(THIS_FOLDER, "ffw")
                run(perf_command + " %s %d %d %d %d" % (executable_abs_path, int(n), fwi, int(l1), int(b)), shell=True)
            else:
                file_abs_path = os.path.join(THIS_FOLDER, fw[fwi]+"_generated_vectorized_tiled.c")
                run("gcc -o ffw-gen " + file_abs_path + flags, shell=True)
                executable_abs_path = os.path.join(THIS_FOLDER, "ffw-gen")
                run(perf_command + " %s %d %d %d" % (executable_abs_path, int(n), int(l1), int(b)), shell=True)


read_config_file()
