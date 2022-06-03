from subprocess import run

import re
import os
from platform import system

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

x = []
y = []
MAX_N = 2048
n = 2048
L1 = 1
B = 1
res = []
tmp_cycles = 0.0
while n <= MAX_N:
    L1 = 1
    while L1 <= n:
        B = 1
        while B <= L1:
            for i in range(repetitions_for_confidence):

                output = run("%s %d %d %d" % (executable_abs_path, n, L1, B), capture_output=True, shell=True).stdout.decode("utf-8")

                cycles = [float(t) for t in output.split() if re.match(r'^-?\d+(?:\.\d+)$', t) is not None][0]
                tmp_cycles = tmp_cycles + cycles

            #Average result over repetitions
            tmp_cycles = tmp_cycles / repetitions_for_confidence
            res.append(((2*(n**3) / tmp_cycles), tmp_cycles, n, L1, B))

            # Print benchmark
            print("Done benchmarking " + str(fw[0]) + " with N = %d, L1 = %d, and B = %d" % (n, L1, B))

            #Update B
            B = B * 2
        #Update L1
        L1 = L1 * 2
    #Update N
    n = n * 2

# Find the best config and write it to a file
with open("result.txt", "w") as file:
    bestRes = []             # Store best results for each n
    cur_n = 4                # Keep track of current n
    best_flops = -1.0          # Keep track of current best perf
    cur_best_config = (0, 0) # L1 & B

    for (flops, cycles, _n, _L1, _B) in res:
        # Find the best config for each n
        if cur_n == _n:
            #print("\nBEST: %lf\nFLOPS: %lf\n" % (best_flops, flops))
            # check for an improved performance
            if flops >= best_flops:
                print("\nOLD_BEST: %lf\nNEW_BEST: %lf\n" % (best_flops, flops))
                best_flops = flops
                cur_best_config = (_L1, _B)
        else:
            # Write the previous best_n to the array
            (best_l1, best_b) = cur_best_config
            bestRes.append((cur_n, best_l1, best_b, best_flops))

            # Prepare next n 
            cur_n = _n
            best_flops = flops
            cur_best_config = (_L1, _B)

        # Write all results to a file and to stdout
        strres = "RESULT: N = %d, L1 = %d, B = %d:\n%d [cycles], %lf [flops/cycle]" % (_n, _L1, _B, cycles, flops)
        print(strres)
        csvres = "%d, %d, %d, %d, %lf\n" % (_n, _L1, _B, cycles, flops)
        file.write(csvres)

    # Write the final best_n to the array
    best_l1, best_b = cur_best_config
    bestRes.append((cur_n, best_l1, best_b, best_flops))

    # Write all best results to a separate file
    with open("best_results.csv", "w") as best_res_file:
        for (best_n, best_l1, best_b, final_flops) in bestRes:
            best_res_file.write("%d, %d, %d, %lf\n" % (best_n, best_l1, best_b, final_flops))
