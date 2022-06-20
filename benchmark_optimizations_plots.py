import os
import string
import sys
from subprocess import run
from typing import Tuple

MAX_N = (2 ** 20)
semi_rings = 3
repetitions_for_confidence = 1

csv_file = "benchmark_results.csv"
fw_base = ["min_plus", "max_min", "or_and"]
implementations = ["baseline", "basic_opt"]
fw = ["min_plus", "or_and", "max_min"]  # functions in generalized floyd warshall
THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
executable_abs_path = os.path.join(THIS_FOLDER, "ffw")


def run_benchmark(fwi: int, n: int, l1: int, b: int) -> Tuple:
    tmp_cycles = 0.0
    for i in range(repetitions_for_confidence):
        output = run("%s %d %d %d %d" % (executable_abs_path, n, fwi, l1, b), capture_output=True,
                     shell=True).stdout.decode("utf-8")

        cycles = float(output.split("\n")[0]) # [float(t) for t in output.split() if re.match(r'^-?\d+(?:\.\d+)$', t) is not None][0]
        tmp_cycles = tmp_cycles + cycles

    # Average result over repetitions
    tmp_cycles = tmp_cycles / repetitions_for_confidence
    res = ((2 * (n ** 3) / tmp_cycles), tmp_cycles, n)

    # Store results in case of crash
    with open(csv_file, "a") as res_dump_file:
        csv_res = "%d,%d,%lf\n" % (n, tmp_cycles, (2 * (n ** 3) / tmp_cycles))
        res_dump_file.write(csv_res)

    # Print benchmark
    print("Done benchmarking " + str(fw[fwi]) + " with N = %d" % (n))
    return res


def set_up(file_name: string):
    file_abs_path = os.path.join(THIS_FOLDER, file_name)
    print(os.getcwd())
    run("gcc -o ffw " + file_abs_path + " tsc_x86.h -march=native -O3 -ffast-math", shell=True)

    for fwi in range(semi_rings):
        print("Benchmarking " + str(fw[fwi]) + " from " + file_name)
        with open(csv_file, "a") as res_dump_file:
            res_dump_file.write(str(fw[fwi]) + " from " + file_name + "\n")

        n = 4096
        res = []

        # Check for config file
        if len(sys.argv) == 2:
            config_file: str = sys.argv[1]
            with open(config_file, "r") as config:
                lines = config.readlines()
                for line in lines:
                    (n, l1, b, _) = line.split(",")
                    print("Begin benchmarking " + str(fw[fwi]) + " from " + file_name + " with N = %s, L1 = %s, and B = "
                                                                                       "%s" % (n, l1, b))
                    res.append(run_benchmark(fwi, int(n), int(l1), int(b)))

        else:
            while n <= MAX_N:
                l1 = n / 8 if (n / 8 > 0) else 1
                b = l1 / 8 if (l1 / 8 > 0) else 1
                res.append(run_benchmark(fwi, n, l1, b))
                # Update N
                n = n * 2

def benchmark_baseline_intermediate(file_name: string):
    file_abs_path = os.path.join(THIS_FOLDER, file_name)
    print(os.getcwd())
    run("gcc -o ffw " + file_abs_path + " tsc_x86.h -march=native -O3 -ffast-math", shell=True)

    for impl in range(2):
        for fwi in range(semi_rings):

            print("Benchmarking " + str(fw_base[fwi]))
            with open(csv_file, "a") as res_dump_file:
                res_dump_file.write(str(fw_base[fwi]) + " " + str(implementations[impl]) + "\n")

            n = 4

            while n <= MAX_N:
                tmp_cycles = 0.0
                for i in range(repetitions_for_confidence):
                    output = run("%s %d %d %d" % (executable_abs_path, n, fwi, impl), capture_output=True, shell=True).stdout.decode("utf-8")
                    cycles = float(output.split()[0])
                    tmp_cycles = tmp_cycles + cycles

                # Average result over repetitions
                tmp_cycles = tmp_cycles / repetitions_for_confidence

                # Store results in case of crash
                with open(csv_file, "a") as res_dump_file:
                    csv_res = "%d, %d, %lf\n" % (n, tmp_cycles, (2*(n**3) / tmp_cycles))
                    res_dump_file.write(csv_res)

                # Print benchmark
                print("Done benchmarking " + str(fw_base[fwi]) + " with N = %d" % (n))

                # Update N
                n = n * 2


# Benchmark baseline and basic optimization
benchmark_baseline_intermediate("basic_optimizations_floyd.c")

# Benchmark tiled version
set_up("tiled_floyd.c")

# Benchmark tiled-vectorized version
set_up("vectorized_tiled_floyd.c")
