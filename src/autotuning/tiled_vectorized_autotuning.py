import os
from subprocess import run
from platform import system

compiled_file = "ffw" + (".exe " if system() == "Windows" else ".out")  # do not remove space

repetitions_for_confidence = 1  # number of repetitions of each n, median is taken
file_name = "../vectorized_tiled_floyd.c"

THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
file_abs_path = os.path.join(THIS_FOLDER, file_name)
executable_abs_path = os.path.join(THIS_FOLDER, "ffw")

fw = ["min_plus", "max_min", "or_and"]  # functions in generalized floyd warshall

run("gcc -o ffw " + file_abs_path + " ../tsc_x86.h -march=native -O3 -ffast-math", shell=True)

MIN_N = 8
MAX_N = (2**12)
MIN_L1 = 4
MAX_L1 = 400
MIN_B = 2
MAX_B = 400

for fwi in range(3):
    with open("autotuned_parameters.csv", "a") as best_res_dump_file:
        best_res_dump_file.write(fw[fwi] + "\n")

    print("Benchmarking " + str(fw[fwi]))
    tmp_cycles = 0.0
    res = []
    n = MIN_N
    while n <= MAX_N:
        res_for_n = []
        L1 = MIN_L1
        while L1 <= n and L1 <= MAX_L1:
            if n % L1 == 0:
                B = MIN_B
                while B <= L1 and B <= MAX_B:
                    if L1 % B == 0:
                        try:
                            tmp_cycles = 0.0
                            for i in range(repetitions_for_confidence):

                                output = run("%s %d %d %d %d" % (executable_abs_path, n, fwi, L1, B),
                                             capture_output=True, shell=True).stdout.decode("utf-8")

                                cycles = float(output.split()[0])
                                tmp_cycles = tmp_cycles + cycles

                            # Average result over repetitions
                            tmp_cycles = tmp_cycles / repetitions_for_confidence
                            res.append(((2*(n**3) / tmp_cycles), tmp_cycles, n, L1, B))
                            res_for_n.append(((2*(n**3) / tmp_cycles), tmp_cycles, n, L1, B))

                            # Store results in case of crash
                            with open("result_dump.csv", "a") as res_dump_file:
                                csv_res = "%s, %d, %d, %d, %d, %lf\n" % (fw[fwi], n, L1, B, tmp_cycles,
                                                                         (2*(n**3) / tmp_cycles))
                                res_dump_file.write(csv_res)

                            # Print benchmark
                            print("Done benchmarking " + str(fw[fwi])
                                  + " with N = %d, L1 = %d, and B = %d" % (n, L1, B))
                        except IndexError:
                            B = B + 1
                            continue

                    # Update B
                    B = B + 1
            # Update L1
            L1 = L1 + 1

        # Compute the best parameters for the past n
        best_flops_for_n = -1.0
        best_config_for_n = (0, 0)
        for (tmp_flops, tmp_cycles, _, tmp_l1, tmp_b) in res_for_n:
            if tmp_flops >= best_flops_for_n:
                best_flops_for_n = tmp_flops
                best_config_for_n = (tmp_l1, tmp_b)

        # Write the best results to a file in case of crash
        with open("autotuned_parameters.csv", "a") as best_res_dump_file:
            (best_l1_for_n, best_b_for_n) = best_config_for_n
            best_res_dump_file.write("%d, %d, %d, %lf\n" % (n, best_l1_for_n, best_b_for_n, best_flops_for_n))

        # Update N
        n = n * 3

    # Compile all of the best configs and write them to a file
    bestRes = []             # Store best results for each n
    cur_n = 4                # Keep track of current n
    best_flops = -1.0        # Keep track of current best perf
    cur_best_config = (0, 0) # L1 & B

    for (flops, cycles, _n, _L1, _B) in res:
        # Find the best config for each n
        if cur_n == _n:
            # check for an improved performance
            if flops >= best_flops:
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
        str_res = "RESULT: N = %d, L1 = %d, B = %d:\n%d [cycles], %lf [flops/cycle]\n" % (_n, _L1, _B, cycles, flops)
        print(str_res)

    # Write the final best_n to the array
    best_l1, best_b = cur_best_config
    bestRes.append((cur_n, best_l1, best_b, best_flops))
