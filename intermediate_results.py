from subprocess import run

import os
from platform import system

compiled_file = "ffw" + (".exe " if system() == "Windows" else ".out") #do not remove space

repetitions_for_confidence = 5 #number of repetitions of each n, median is taken
semi_rings = 3
file_name = "basic_optimizations_floyd.c" #do not remove space after the end!

THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
file_abs_path = os.path.join(THIS_FOLDER, file_name)
executable_abs_path = os.path.join(THIS_FOLDER, "ffw")

print(os.getcwd())
fw = ["min_plus", "max_min", "or_and"] #functions in generalized floyd warshall

run("gcc -o ffw " + file_abs_path + " tsc_x86.h -O0", shell=True)
# run("gcc -o ffw " + file_abs_path + " tsc_x86.h -march=native -O3 -ffast-math", shell=True)


for j in range(semi_rings):

    print("Benchmarking " + str(fw[j]))
    with open("flagO0_basic_opt.csv", "a") as res_dump_file:
        res_dump_file.write(str(fw[j]) + "\n")

    MAX_N = (256)
    n = 4
    tmp_cycles = 0.0

    while n <= MAX_N:
        tmp_cycles = 0.0
        for i in range(repetitions_for_confidence):

                output = run("%s %d %d" % (executable_abs_path, n, j), capture_output=True, shell=True).stdout.decode("utf-8")

                cycles = float(output.split()[0])
                tmp_cycles = tmp_cycles + cycles

        #Average result over repetitions
        tmp_cycles = tmp_cycles / repetitions_for_confidence

        # Store results in case of crash
        with open("flagO0_basic_opt.csv", "a") as res_dump_file:
            csvres = "%d, %d, %lf\n" % (n, tmp_cycles, (2*(n**3) / tmp_cycles))
            res_dump_file.write(csvres)

        # Print benchmark
        print("Done benchmarking " + str(fw[j]) + " with N = %d" % (n))

        #Update N
        n = n * 2
