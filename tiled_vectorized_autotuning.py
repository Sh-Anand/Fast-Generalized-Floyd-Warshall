from subprocess import run

import re
import os
from platform import system

compiled_file = "ffw" + (".exe " if system() == "Windows" else ".out") #do not remove space

repetitions_for_confidence = 10 #number of repetitions of each n, median is taken
file_name = "vectorized_tiled_floyd.c" #do not remove space after the end!

THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
file_abs_path = os.path.join(THIS_FOLDER, file_name)
executable_abs_path = os.path.join(THIS_FOLDER, "ffw")

print(os.getcwd())
fw = ["min_plus", "or_and", "max_min"] #functions in generalized floyd warshall

print("Benchmarking " + str(fw[0]))

run("gcc -o ffw " + file_abs_path + " tsc_x86.h -march=native -O3 -ffast-math", shell=True)

x = []
y = []
MAX_N = 1024
n = 4
L1 = 1
B = 1
res = []
tmp_cycles = 0.0
while n <= MAX_N:
    while L1 <= n:
        while B <= L1:
            for i in range(repetitions_for_confidence):

                output = run("%s %d %d %d" % (executable_abs_path, n, L1, B), capture_output=True).stdout.decode("utf-8")

                cycles = [float(t) for t in output.split() if re.match(r'^-?\d+(?:\.\d+)$', t) is not None][0]
                tmp_cycles = tmp_cycles + cycles

            #Average result over repetitions
            tmp_cycles = tmp_cycles / repetitions_for_confidence
            res.append(((2*(n**3) / tmp_cycles), tmp_cycles, n, L1, B))

            #Update B
            B = B * 2
        #Update L1
        L1 = L1 * 2
    #Update N
    n = n * 2

    print("Done benchmarking " + str(fw[0]))

for (flops, cycles, _n, _L1, _B) in res:
    print("RESULT: N = %d, L1 = %d, B = %d:\n%d [cycles], %lf [flops/cycle]" % (_n, _L1, _B, cycles, flops))