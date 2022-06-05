from subprocess import run

import re
import os
from platform import system
import sys
from typing import Tuple

from pytest import param

def run_benchmark(n: int, l1: int, b: int) -> Tuple:
    repetitions_for_confidence = 5
    tmp_cycles = 0.0
    for i in range(repetitions_for_confidence):

        output = run("%s %d %d %d" % (executable_abs_path, n, l1, b), capture_output=True, shell=True).stdout.decode("utf-8")

        cycles = float(output.split("\n")[0])#[float(t) for t in output.split() if re.match(r'^-?\d+(?:\.\d+)$', t) is not None][0]
        tmp_cycles = tmp_cycles + cycles

    #Average result over repetitions
    tmp_cycles = tmp_cycles / repetitions_for_confidence
    res = ((2*(n**3) / tmp_cycles), tmp_cycles, n)

    # Store results in case of crash
    with open("result_dump_no_opt.csv", "a") as res_dump_file:
        csvres = "%d,%d,%lf\n" % (n, tmp_cycles, (2*(n**3) / tmp_cycles))
        res_dump_file.write(csvres)

    # Print benchmark
    print("Done benchmarking " + str(fw[0]) + " with N = %d" % (n))
    return res

compiled_file = "ffw" + (".exe " if system() == "Windows" else ".out") #do not remove space

repetitions_for_confidence = 1 #number of repetitions of each n, median is taken
file_name = "vectorized_tiled_floyd.c" #do not remove space after the end!

THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
file_abs_path = os.path.join(THIS_FOLDER, file_name)
executable_abs_path = os.path.join(THIS_FOLDER, "ffw")

print(os.getcwd())
fw = ["min_plus", "or_and", "max_min"] #functions in generalized floyd warshall

print("Benchmarking " + str(fw[0]))

run("gcc -o ffw " + file_abs_path + " tsc_x86.h -march=native -O3 -ffast-math", shell=True)

MAX_N = (2**20)
n = 4096
res = []
tmp_cycles = 0.0

# Check for config file
if len(sys.argv) == 2:
    config_file: str = sys.argv[1]
    with open(config_file, "r") as config:
        lines = config.readlines()
        for line in lines:
            (n, l1, b, _) = line.split(",")
            print("Begin benchmarking " + str(fw[0]) + " with N = %s, L1 = %s, and B = %s" % (n, l1, b))
            res.append(run_benchmark(int(n), int(l1), int(b)))

else:
    while n <= MAX_N:
        l1 = n/8 if (n/8 > 0) else 1
        b = l1/8 if (l1/8 > 0) else 1
        res.append(run_benchmark(n, l1, b))

        #Update N
        n = n * 2


    

