import os
import string
import sys
from subprocess import run

MIN_N = 8
MAX_N = (2 ** 12)
semi_rings = 3
repetitions_for_confidence = 10

csv_file = "benchmark_results.csv"
fw = ["min_plus", "max_min", "or_and"]  # functions in generalized floyd warshall
implementations = ["baseline", "basic_opt"]
THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
executable_abs_path = os.path.join(THIS_FOLDER, "ffw")


def run_benchmark(generated: bool, fwi: int, n: int, l1: int, b: int):
    tmp_cycles = 0.0
    for i in range(repetitions_for_confidence):
        if generated:
            if n < 216:     # If small N, run without unrolling
                file_name = "../vectorized_tiled_floyd.c"
                file_abs_path = os.path.join(THIS_FOLDER, file_name)
                run("gcc -o ffw1 " + file_abs_path + " ../tsc_x86.h -march=native -O3 -ffast-math", shell=True)
                output = run("%s %d %d %d %d" % (os.path.join(THIS_FOLDER, "ffw1"), n, fwi, l1, b), capture_output=True,
                             shell=True).stdout.decode("utf-8")
            else:
                output = run("%s %d %d %d" % (executable_abs_path, n, l1, b), capture_output=True,
                         shell=True).stdout.decode("utf-8")
        else:
            output = run("%s %d %d %d %d" % (executable_abs_path, n, fwi, l1, b), capture_output=True,
                         shell=True).stdout.decode("utf-8")

        cycles = float(output.split("\n")[0])
        tmp_cycles = tmp_cycles + cycles

    # Average result over repetitions
    tmp_cycles = tmp_cycles / repetitions_for_confidence

    # Store results in case of crash
    with open(csv_file, "a") as res_dump_file:
        csv_res = "%d,%d,%lf\n" % (n, tmp_cycles, (2 * (n ** 3) / tmp_cycles))
        res_dump_file.write(csv_res)

    # Print benchmark
    print("Done benchmarking " + str(fw[fwi]) + " with N = %d" % n)


def set_up(file_name: string, generated: bool):
    # Check for config file
    if len(sys.argv) != 2:
        print("Error: Configuration file needed")
        return

    if not generated:
        file_abs_path = os.path.join(THIS_FOLDER, file_name)
        run("gcc -o ffw " + file_abs_path + " ../tsc_x86.h -march=native -O3 -ffast-math", shell=True)
    fwi = -1

    config_file: str = sys.argv[1]
    with open(config_file, "r") as config:
        lines = config.readlines()
        for line in lines:
            # Handle each semi-ring
            if len(line.split(",")) < 2:
                fwi = fwi + 1
                # If benchmarking generated code, run corresponding file
                if generated:
                    file_name = fw[fwi]+"_generated_vectorized_tiled.c"
                    file_abs_path = os.path.join(THIS_FOLDER, file_name)
                    run("gcc -o ffw " + file_abs_path + " ../tsc_x86.h -march=native -O3 -ffast-math", shell=True)

                with open(csv_file, "a") as res_dump_file:
                    res_dump_file.write(str(fw[fwi]) + " " + file_name + "\n")

                continue

            (n, l1, b, _) = line.split(",")
            print("Begin benchmarking " + str(fw[fwi]) + " from " + file_name +
                  " with N = %s, L1 = %s, and B = %s" % (n, l1, b))
            run_benchmark(generated, fwi, int(n), int(l1), int(b))


def benchmark_baseline_intermediate(file_name: string):
    file_abs_path = os.path.join(THIS_FOLDER, file_name)
    run("gcc -o ffw " + file_abs_path + " ../tsc_x86.h -march=native -O3 -ffast-math", shell=True)

    for impl in range(2):
        for fwi in range(semi_rings):

            print("Benchmarking " + str(fw[fwi]))
            with open(csv_file, "a") as res_dump_file:
                res_dump_file.write(str(fw[fwi]) + " " + str(implementations[impl]) + "\n")

            n = MIN_N

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
                print("Done benchmarking " + str(fw[fwi]) + " with N = %d" % n)

                # Update N
                n = n * 3


def run_generator():
    for fwi in range(semi_rings):
        # Generate unrolled code for each semi-ring
        run("python3 ../generator/generator.py " + str(fwi), shell=True)


# Benchmark baseline and basic optimization
benchmark_baseline_intermediate("../basic_optimizations_floyd.c")

# Benchmark tiled version
set_up("../tiled_floyd.c", False)

# Benchmark tiled-vectorized version
set_up("../vectorized_tiled_floyd.c", False)

# Benchmark generated unrolled tiled-vectorized version
run_generator()
set_up("", True)
